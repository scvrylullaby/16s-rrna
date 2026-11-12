import type { Result } from './ncbi.ts'
import { ambiguous, emptyResult, hash, validSequence } from './ncbi.ts'
import { strains } from './strains.ts'

const columns = [
  'id',
  'name',
  'status',
  'accession',
  'organism',
  'sourceNames',
  'source',
  'location',
  'length',
  'ambiguousBases',
  'rrnaStatus',
  'checkedAt',
  'note',
  'query',
  'url',
  'sha256'
] as const
export interface Dataset {
  results: Result[]
  sequences: Map<string, string>
}
const file = (directory: string | URL, name: string): string | URL =>
  directory instanceof URL
    ? new URL(name, directory)
    : `${directory.replace(/[\\/]+$/, '')}/${name}`

// CSV допускает запятые, кавычки и переносы строк внутри заключённого в кавычки поля.
export function parseCsv(csv: string): string[][] {
  const rows: string[][] = []
  let row: string[] = [],
    field = '',
    quoted = false
  csv = csv.replace(/^\uFEFF/, '')
  for (let index = 0; index < csv.length; index++) {
    const char = csv[index]
    if (char === '"') {
      if (quoted && csv[index + 1] === '"') {
        field += '"'
        index++
      } else quoted = !quoted
    } else if (char === ',' && !quoted) {
      row.push(field)
      field = ''
    } else if ((char === '\r' || char === '\n') && !quoted) {
      if (char === '\r' && csv[index + 1] === '\n') index++
      row.push(field)
      if (row.some(value => value !== '')) rows.push(row)
      row = []
      field = ''
    } else field += char
  }
  if (quoted) throw new Error('Незакрытые кавычки в CSV')
  row.push(field)
  if (row.some(value => value !== '')) rows.push(row)
  return rows
}

export function parseFasta(fasta: string): Map<string, string> {
  const sequences = new Map<string, string>()
  let id = ''
  for (const line of fasta.split(/\r?\n/).filter(Boolean)) {
    if (line.startsWith('>')) {
      id = line.slice(1).split('|')[0]
      if (!/^\d{2}$/.test(id) || sequences.has(id)) {
        throw new Error('Некорректный или повторяющийся ID в FASTA')
      }
      sequences.set(id, '')
    } else {
      if (!id) throw new Error('Последовательность без заголовка FASTA')
      sequences.set(id, sequences.get(id)! + line.trim())
    }
  }
  return sequences
}

export async function validateDataset(dataset: Dataset): Promise<number> {
  if (dataset.results.length !== strains.length) {
    throw new Error('В таблице должно быть 30 штаммов')
  }
  const ids = new Set<string>()
  for (const result of dataset.results) {
    if (ids.has(result.id)) throw new Error(`Повтор ID ${result.id}`)
    ids.add(result.id)
    const strain = strains.find(item => item.id === result.id)
    if (!strain || strain.name !== result.name)
      throw new Error(`Неверное имя штамма ${result.id}`)
    const sequence = dataset.sequences.get(result.id)
    if (result.status !== 'downloaded') {
      if (sequence !== undefined)
        throw new Error(`Лишняя последовательность ${result.id}`)
      continue
    }
    if (
      !sequence ||
      !validSequence(sequence) ||
      sequence.length !== result.length ||
      ambiguous(sequence) !== result.ambiguousBases ||
      (await hash(sequence)) !== result.sha256
    ) {
      throw new Error(`FASTA не соответствует таблице: ${result.id}`)
    }
    if (
      !result.accession ||
      !['rrna_record', 'gene_record', 'genome_feature'].includes(result.source)
    ) {
      throw new Error(`Нет источника последовательности ${result.id}`)
    }
    if (
      strain.identityAccessions.length &&
      !strain.identityAccessions.some(
        accession =>
          accession.replace(/^NZ_/, '').split('.')[0] ===
          result.accession.replace(/^NZ_/, '').split('.')[0]
      )
    ) {
      throw new Error(
        `Не подтверждена культура ${result.id}: ${result.accession}`
      )
    }
  }
  if ([...dataset.sequences.keys()].some(id => !ids.has(id))) {
    throw new Error('Неизвестный штамм в FASTA')
  }
  return dataset.sequences.size
}

export async function loadDataset(directory: string | URL): Promise<Dataset> {
  let csv: string
  try {
    csv = await Deno.readTextFile(file(directory, 'results.csv'))
  } catch (error) {
    if (!(error instanceof Deno.errors.NotFound)) throw error
    return { results: strains.map(emptyResult), sequences: new Map() }
  }
  const [header, ...rows] = parseCsv(csv)
  if (JSON.stringify(header) !== JSON.stringify(columns)) {
    throw new Error('Неожиданный формат таблицы results.csv')
  }
  const results = rows.map(values => {
    if (values.length !== columns.length)
      throw new Error('Некорректная строка CSV')
    const row = Object.fromEntries(
      columns.map((key, index) => [key, values[index]])
    ) as Record<keyof Result, string>
    if (
      !['downloaded', 'not_found', 'needs_review', 'access_error'].includes(
        row.status
      )
    ) {
      throw new Error('Неизвестный статус')
    }
    return {
      ...row,
      status: row.status as Result['status'],
      source: row.source as Result['source'],
      length: row.length === '' ? null : Number(row.length),
      ambiguousBases:
        row.ambiguousBases === '' ? null : Number(row.ambiguousBases)
    }
  })
  const sequences = parseFasta(
    await Deno.readTextFile(file(directory, '16s_all.fasta'))
  )
  const dataset = { results, sequences }
  await validateDataset(dataset)
  const strict = parseFasta(
    await Deno.readTextFile(file(directory, '16s_rrna.fasta'))
  )
  const expected = results.filter(
    row => row.status === 'downloaded' && row.source === 'rrna_record'
  )
  if (
    strict.size !== expected.length ||
    expected.some(row => strict.get(row.id) !== sequences.get(row.id))
  ) {
    throw new Error('FASTA с фильтром rRNA не соответствует таблице')
  }
  return dataset
}

export async function saveDataset(
  directory: string | URL,
  dataset: Dataset
): Promise<void> {
  await validateDataset(dataset)
  await Deno.mkdir(directory, { recursive: true })
  const escape = (value: string | number | null) => {
    const text = String(value ?? '')
    return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text
  }
  const csv =
    '\uFEFF' +
    [
      columns.join(','),
      ...dataset.results.map(result =>
        columns.map(key => escape(result[key])).join(',')
      )
    ].join('\n') +
    '\n'
  const all: string[] = [],
    strict: string[] = []
  for (const result of dataset.results) {
    const sequence = dataset.sequences.get(result.id)
    if (!sequence) continue
    const header = `>${result.id}|${result.accession}|${result.name.replaceAll(
      ' ',
      '_'
    )}|16S|location=${result.location}|source=${result.source}`
    const fasta = header + '\n' + sequence.match(/.{1,70}/g)!.join('\n') + '\n'
    all.push(fasta)
    if (result.source === 'rrna_record') strict.push(fasta)
  }
  for (const [filename, content] of [
    ['results.csv', csv],
    ['16s_all.fasta', all.join('')],
    ['16s_rrna.fasta', strict.join('')]
  ]) {
    await Deno.writeTextFile(file(directory, filename + '.tmp'), content)
    await Deno.rename(
      file(directory, filename + '.tmp'),
      file(directory, filename)
    )
  }
}

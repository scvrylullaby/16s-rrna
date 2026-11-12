import { DOMParser } from '@xmldom/xmldom'
import type { Strain } from './strains.ts'

export type Source = 'rrna_record' | 'gene_record' | 'genome_feature' | ''
export type Status =
  | 'downloaded'
  | 'not_found'
  | 'needs_review'
  | 'access_error'
export interface Result {
  id: string
  name: string
  status: Status
  accession: string
  organism: string
  sourceNames: string
  source: Source
  location: string
  length: number | null
  ambiguousBases: number | null
  rrnaStatus: string
  checkedAt: string
  note: string
  query: string
  url: string
  sha256: string
}
export interface Candidate {
  accession: string
  organism: string
  sourceNames: string
  location: string
  sequence: string
  reason: string
  range?: { start: number; stop: number; reverse: boolean }
}
export type XmlRequest = (
  endpoint: string,
  params: Record<string, string | number>
) => Promise<Document>

export const elements = (node: Document | Element, tag: string): Element[] =>
  Array.from(node.getElementsByTagName(tag))
export const text = (node: Document | Element, tag: string): string =>
  elements(node, tag)[0]?.textContent ?? ''
const normalize = (value: string) =>
  value.toLowerCase().replace(/[^a-z0-9]/g, '')
const baseAccession = (value: string) => value.replace(/^NZ_/, '').split('.')[0]
export async function hash(sequence: string): Promise<string> {
  const digest = await crypto.subtle.digest(
    'SHA-256',
    new TextEncoder().encode(sequence)
  )
  return Array.from(new Uint8Array(digest), byte =>
    byte.toString(16).padStart(2, '0')
  ).join('')
}
export const ambiguous = (sequence: string) =>
  (sequence.match(/[^ACGT]/g) ?? []).length
export const validSequence = (sequence: string) =>
  /^[ACGTRYMKSWBDHVN]+$/.test(sequence)

export function parseXml(xml: string): Document {
  const fail = (message: unknown): never => {
    throw new Error(`Некорректный XML: ${message}`)
  }
  const document = new DOMParser({
    errorHandler: { warning: fail, error: fail, fatalError: fail }
  }).parseFromString(xml, 'application/xml')
  if (!document.documentElement) throw new Error('Пустой XML')
  const errors = [
    ...elements(document, 'ERROR'),
    ...elements(document, 'ErrorList')
  ]
  if (errors.length)
    throw new Error(errors.map(node => node.textContent).join('; '))
  return document
}

export function reverseComplement(sequence: string): string {
  const pairs: Record<string, string> = {
    A: 'T',
    C: 'G',
    G: 'C',
    T: 'A',
    R: 'Y',
    Y: 'R',
    M: 'K',
    K: 'M',
    S: 'S',
    W: 'W',
    B: 'V',
    V: 'B',
    D: 'H',
    H: 'D',
    N: 'N'
  }
  return [...sequence]
    .reverse()
    .map(base => pairs[base])
    .join('')
}

// Проверяем культуру по source, затем извлекаем именно аннотированный ген 16S.
export function checkRecord(record: Element, strain: Strain): Candidate {
  const result: Candidate = {
    accession: text(record, 'GBSeq_accession-version'),
    organism: text(record, 'GBSeq_organism'),
    sourceNames: '',
    location: '',
    sequence: '',
    reason: ''
  }
  const reject = (reason: string) => ({ ...result, reason })
  if (!/^(Rhodococcus|Prescottella) /.test(result.organism))
    return reject('Другой род организма')
  if (
    strain.identityAccessions.length &&
    !strain.identityAccessions.some(
      value => baseAccession(value) === baseAccession(result.accession)
    )
  ) {
    return reject(
      'Короткое имя штамма неоднозначно; accession не связан с нужной культурой'
    )
  }
  if (
    /record[\s\S]{0,30}(suppressed|withdrawn|removed)/i.test(
      text(record, 'GBSeq_comment')
    )
  ) {
    return reject('Запись отозвана')
  }
  const features = elements(record, 'GBFeature')
  const sourceNames = features
    .filter(feature => text(feature, 'GBFeature_key') === 'source')
    .flatMap(feature => elements(feature, 'GBQualifier'))
    .filter(qualifier =>
      ['strain', 'isolate', 'culture_collection'].includes(
        text(qualifier, 'GBQualifier_name')
      )
    )
    .flatMap(qualifier => text(qualifier, 'GBQualifier_value').split(/[;=()]/))
    .map(name => name.trim())
    .filter(Boolean)
  result.sourceNames = sourceNames.join('; ')
  if (
    !sourceNames.some(name =>
      strain.aliases.some(alias => normalize(alias) === normalize(name))
    )
  ) {
    return reject('Обозначение штамма не подтверждено в source')
  }
  const sequence = text(record, 'GBSeq_sequence').toUpperCase()
  const candidates: Candidate[] = []
  for (const feature of features) {
    if (text(feature, 'GBFeature_key') !== 'rRNA') continue
    const products = elements(feature, 'GBQualifier')
      .filter(qualifier => text(qualifier, 'GBQualifier_name') === 'product')
      .map(qualifier => text(qualifier, 'GBQualifier_value'))
    if (
      !products.some(
        product =>
          /\b16S\b/i.test(product) && !/spacer|intergenic/i.test(product)
      )
    )
      continue
    const intervals = elements(feature, 'GBInterval')
    if (intervals.length !== 1) continue // Составные координаты требуют ручной проверки.
    const interval = intervals[0]
    if (
      (text(interval, 'GBInterval_accession') || result.accession) !==
      result.accession
    )
      continue
    const from = Number(text(interval, 'GBInterval_from'))
    const to = Number(text(interval, 'GBInterval_to'))
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      Math.min(from, to) < 1
    )
      continue
    const start = Math.min(from, to),
      stop = Math.max(from, to)
    const reverse =
      from > to ||
      elements(interval, 'GBInterval_iscomp')[0]?.getAttribute('value') ===
        'true'
    const candidate = {
      ...result,
      location: text(feature, 'GBFeature_location')
    }
    if (!sequence) {
      candidates.push({ ...candidate, range: { start, stop, reverse } })
    } else if (stop <= sequence.length) {
      let fragment = sequence.slice(start - 1, stop) // NCBI использует координаты с 1, обе границы включены.
      if (reverse) fragment = reverseComplement(fragment)
      if (validSequence(fragment))
        candidates.push({ ...candidate, sequence: fragment })
    }
  }
  const length = (candidate: Candidate) =>
    candidate.sequence.length ||
    (candidate.range ? candidate.range.stop - candidate.range.start + 1 : 0)
  candidates.sort(
    (a, b) =>
      length(b) - length(a) || ambiguous(a.sequence) - ambiguous(b.sequence)
  )
  return (
    candidates[0] ?? reject('Нет доступного аннотированного фрагмента 16S rRNA')
  )
}

export function buildQueries(
  strain: Strain
): { source: Source; query: string }[] {
  const terms = strain.aliases
    .map(alias => `"${alias}"[All Fields]`)
    .join(' OR ')
  const taxa = strain.taxa.map(taxon => `"${taxon}"[Organism]`).join(' OR ')
  const base = `(${taxa}) AND (${terms})`
  const known = strain.accessions
    .map(value => `"${value}"[Accession]`)
    .join(' OR ')
  const anchors = strain.identityAccessions
    .map(value => `"${value}"[Accession]`)
    .join(' OR ')
  return [
    {
      source: 'rrna_record',
      query: `${base} AND biomol_rrna[PROP] AND 16S[Title]`
    },
    {
      source: 'gene_record',
      query: `${known ? `(${base} OR (${known}))` : base} AND 16S[Title]`
    },
    {
      source: 'genome_feature',
      query: `${
        anchors ? `(${anchors})` : base
      } AND biomol_genomic[PROP] AND rrna[FKEY] NOT plasmid[Title]`
    }
  ]
}

export async function requestXml(
  endpoint: string,
  params: Record<string, string | number>
): Promise<Document> {
  const query = new URLSearchParams(
    Object.entries({ ...params, tool: 'student_16s_lab' }).map(
      ([key, value]) => [key, String(value)]
    )
  )
  const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/${endpoint}?${query}`
  for (let attempt = 0; attempt < 3; attempt++) {
    await new Promise(resolve =>
      setTimeout(resolve, attempt === 0 ? 400 : 2000 * attempt)
    )
    let response: Response
    try {
      response = await fetch(url, {
        signal: AbortSignal.timeout(60_000),
        headers: { 'User-Agent': 'student-16s-lab/2.0' }
      })
    } catch (error) {
      if (attempt === 2) throw error
      continue
    }
    if (response.ok) return parseXml(await response.text())
    if (![429, 500, 502, 503, 504].includes(response.status) || attempt === 2) {
      throw new Error(`NCBI HTTP ${response.status}`)
    }
  }
  throw new Error('NCBI недоступен')
}

export function emptyResult(strain: Strain): Result {
  return {
    id: strain.id,
    name: strain.name,
    status: 'needs_review',
    accession: '',
    organism: '',
    sourceNames: '',
    source: '',
    location: '',
    length: null,
    ambiguousBases: null,
    rrnaStatus: 'not_searched',
    checkedAt: '',
    note: 'Поиск ещё не выполнен',
    query: '',
    url: '',
    sha256: ''
  }
}

export async function collect(
  strain: Strain,
  request: XmlRequest = requestXml
): Promise<{ result: Result; sequence: string }> {
  const result = {
    ...emptyResult(strain),
    checkedAt: new Date().toISOString(),
    note: ''
  }
  const seen = new Set<string>()
  let hadCandidates = false,
    hadWarnings = false
  try {
    for (const stage of buildQueries(strain)) {
      result.query = stage.query
      const ids: string[] = [],
        warnings: string[] = []
      let count = 0,
        offset = 0
      do {
        const search = await request('esearch.fcgi', {
          db: 'nuccore',
          term: stage.query,
          retstart: offset,
          retmax: 200
        })
        count = Number(text(search, 'Count'))
        if (
          !Number.isInteger(count) ||
          count < 0 ||
          !elements(search, 'Count').length
        ) {
          throw new Error('Нет корректного Count в ответе NCBI')
        }
        ids.push(...elements(search, 'Id').map(node => node.textContent ?? ''))
        warnings.push(
          ...elements(search, 'WarningList').map(node => node.textContent ?? '')
        )
        offset += 200
      } while (offset < count)
      if (ids.length !== count || new Set(ids).size !== count) {
        throw new Error('Выдача NCBI получена не полностью')
      }
      hadCandidates ||= count > 0
      hadWarnings ||= warnings.length > 0
      if (stage.source === 'rrna_record') {
        result.rrnaStatus = count
          ? 'no_exact_match'
          : warnings.length
            ? 'query_warning'
            : 'zero_results'
      }
      const candidates: Candidate[] = []
      for (const id of ids) {
        if (seen.has(id)) continue
        seen.add(id)
        const fetched = await request('efetch.fcgi', {
          db: 'nuccore',
          id,
          rettype: 'gb',
          retmode: 'xml'
        })
        const records = elements(fetched, 'GBSeq')
        if (records.length !== 1)
          throw new Error('NCBI вернул неожиданное количество записей')
        const candidate = checkRecord(records[0], strain)
        if (candidate.reason) continue
        if (candidate.range) {
          const { start, stop, reverse } = candidate.range
          const slice = await request('efetch.fcgi', {
            db: 'nuccore',
            id: candidate.accession,
            rettype: 'gb',
            retmode: 'xml',
            seq_start: start,
            seq_stop: stop,
            strand: reverse ? 2 : 1
          })
          if (text(slice, 'GBSeq_accession-version') !== candidate.accession) {
            throw new Error('Получен другой accession фрагмента')
          }
          candidate.sequence = text(slice, 'GBSeq_sequence').toUpperCase()
          if (
            candidate.sequence.length !== stop - start + 1 ||
            !validSequence(candidate.sequence)
          ) {
            throw new Error(
              'Последовательность не соответствует координатам 16S'
            )
          }
        }
        if (candidate.sequence) candidates.push(candidate)
      }
      candidates.sort(
        (a, b) =>
          b.sequence.length - a.sequence.length ||
          ambiguous(a.sequence) - ambiguous(b.sequence)
      )
      const best = candidates[0]
      if (!best) continue
      Object.assign(result, {
        status: 'downloaded',
        accession: best.accession,
        organism: best.organism,
        sourceNames: best.sourceNames,
        source: stage.source,
        location: best.location,
        length: best.sequence.length,
        ambiguousBases: ambiguous(best.sequence),
        url: `https://www.ncbi.nlm.nih.gov/nuccore/${best.accession}`,
        sha256: await hash(best.sequence),
        note:
          stage.source === 'rrna_record' ? '' : 'Дополнительный поиск по ДНК'
      })
      if (stage.source === 'rrna_record') result.rrnaStatus = 'downloaded'
      return { result, sequence: best.sequence }
    }
    result.status = hadCandidates || hadWarnings ? 'needs_review' : 'not_found'
    result.note =
      hadCandidates || hadWarnings
        ? 'Кандидаты или предупреждения требуют проверки'
        : 'Запросы дали 0 записей; это не доказательство отсутствия гена'
  } catch (error) {
    result.status = 'access_error'
    result.note = error instanceof Error ? error.message : String(error)
    if (result.rrnaStatus === 'not_searched') result.rrnaStatus = 'access_error'
  }
  return { result, sequence: '' }
}

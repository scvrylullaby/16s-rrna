import { collect } from './ncbi.ts';
import { loadDataset, saveDataset } from './data.ts';
import { strains } from './strains.ts';

try {
  const directory = new URL('../data/', import.meta.url);
  const dataset = await loadDataset(directory);
  for (const strain of strains) {
    console.log(`${strain.id}. ${strain.name}: поиск…`);
    const { result, sequence } = await collect(strain);
    const index = dataset.results.findIndex((row) => row.id === strain.id);
    if (result.status === 'downloaded') {
      dataset.results[index] = result;
      dataset.sequences.set(strain.id, sequence);
    } else if (dataset.results[index].status !== 'downloaded') {
      dataset.results[index] = result;
    } else console.log('  Сохраняю ранее подтверждённую последовательность.');
    await saveDataset(directory, dataset);
    console.log(
      `  ${result.status}${
        result.accession ? `: ${result.accession}, ${result.length} нт` : `: ${result.note}`
      }`,
    );
    if (result.status === 'access_error' || result.status === 'needs_review') Deno.exitCode = 1;
    if (result.status === 'access_error') break;
  }
  console.log(`Сохранено: ${dataset.sequences.size}/30. Файлы: data/`);
} catch (error) {
  console.error(error instanceof Error ? error.message : String(error));
  Deno.exitCode = 1;
}

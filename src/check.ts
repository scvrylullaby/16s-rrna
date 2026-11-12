import { loadDataset } from './data.ts';

try {
  const dataset = await loadDataset(new URL('../data/', import.meta.url));
  console.log(`Таблица и FASTA согласованы: ${dataset.sequences.size} последо вательностей из 30.`);
  for (const row of dataset.results.filter((row) => row.status !== 'downloaded')) {
    console.log(`${row.id}. ${row.name}: ${row.note}`);
  }
} catch (error) {
  console.error(error instanceof Error ? error.message : String(error));
  Deno.exitCode = 1;
}

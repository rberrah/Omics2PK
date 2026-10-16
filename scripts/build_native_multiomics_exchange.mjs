#!/usr/bin/env node
/**
 * Generate the REAL MultiOmics V1 exchange with its original JS producer.
 * Inputs are the publicly shipped, SYNTHETIC three-omic demo CSVs only.
 * No patient records, clinical data or network inference are processed.
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, '../upstream/multiomics');
const source = await import(path.join(root, 'src/lib/multiomics/downstream-exchange.js'));
const demo = path.join(root, 'static/multiomics');
function read(name) {
  return fs.readFileSync(path.join(demo, name), 'utf8');
}
const inputs = {
  metadata: { content: read('demo_metadata.csv') },
  transcriptomics: { content: read('demo_transcriptomics.csv') },
  proteomics: { content: read('demo_proteomics.csv') },
  metabolomics: { content: read('demo_metabolomics.csv') },
};
const project = {
  inputs,
  settings: {
    unitType: 'participant',
    organism: 'human',
    transcriptomicsValues: 'raw_counts',
    transcriptomicsIdType: 'unresolved_source_namespace',
    proteomicsValues: 'source_values_unspecified',
    proteomicsIdType: 'uniprot',
    metabolomicsValues: 'source_values_unspecified',
    metabolomicsIdType: 'chebi',
  },
};
const result = await source.buildDownstreamExchange(project);
const output = process.argv[2];
if (!output) throw new Error('Provide a local output directory');
const resolved = path.resolve(output);
fs.mkdirSync(resolved, {recursive: true});
for (const item of result.entries) {
  const rel = item.name.replace(/^exchange\//, '');
  if (rel === item.name || path.isAbsolute(rel) || rel.includes('..')) {
    throw new Error('Unexpected producer output path');
  }
  const dest = path.join(resolved, rel);
  fs.mkdirSync(path.dirname(dest), {recursive: true});
  fs.writeFileSync(dest, item.data, 'utf8');
}
console.log(JSON.stringify({
  source: 'MultiOmics native browser source export',
  demo: 'repository bundled SYNTHETIC demos only',
  files: result.entries.length,
  omics: result.manifest.matrices.map(x => x.omic),
  status: 'BUILT_NOT_YET_VALIDATED',
}));

#!/usr/bin/env node
/** Original MultiOmics research ZIP PACKER, using bundled synthetic CSVs only.
 * This exercises the end-user ZIP layer (not just the nested exchange producer).
 * The synthetic result is explicitly a tiny test stub, NOT a biological fit.
 */
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const here=path.dirname(fileURLToPath(import.meta.url));
const upstream=path.resolve(here,'../upstream/multiomics');
const {buildResearchBundle}=await import(path.join(upstream,'src/lib/multiomics/research-export-package.js'));
const demo=path.join(upstream,'static/multiomics');
const read=name=>fs.readFileSync(path.join(demo,name),'utf8');
const inputs={
  metadata:{content:read('demo_metadata.csv'),sha256:'synthetic-demo-source-metadata'},
  transcriptomics:{content:read('demo_transcriptomics.csv'),sha256:'synthetic-demo-source-transcriptomics'},
  proteomics:{content:read('demo_proteomics.csv'),sha256:'synthetic-demo-source-proteomics'},
  metabolomics:{content:read('demo_metabolomics.csv'),sha256:'synthetic-demo-source-metabolomics'},
};
const project={
  format:'pmx-multiomics-portable-project',
  schemaVersion:1,
  inputs,
  settings:{
    unitType:'participant',organism:'human',
    transcriptomicsValues:'raw_counts',transcriptomicsIdType:'unresolved_source_namespace',
    proteomicsValues:'source_values_unspecified',proteomicsIdType:'uniprot',
    metabolomicsValues:'source_values_unspecified',metabolomicsIdType:'chebi',
  },
  archivedResult:null
};
const result={
  generatedAt:'2026-10-16T00:00:00.000Z',
  engine:{name:'synthetic-research-ZIP-integration-test-only',version:'0'},
  protocol:{objective:'explore',designType:'repeated',analysisIntent:'technical_fixture'},
  layers:{transcriptomics:{rows:[{feature:'SYNTH_SOURCE_FEATURE',effect:0,pValue:null,qValue:null}],inferenceMethod:'fake_test_stub_no_science'}},
};
const args=process.argv.slice(2);
if(args.length!==1)throw Error('Usage node build_native_multiomics_research_zip.mjs /tmp/output.zip');
const bundle=await buildResearchBundle({result,portableProject:project,locale:'fr'});
if(bundle.manifest.downstream_exchange.status!=='included')throw Error('Original package did not contain real V1 exchange: '+JSON.stringify(bundle.manifest.downstream_exchange));
const out=path.resolve(args[0]);
fs.mkdirSync(path.dirname(out),{recursive:true});
fs.writeFileSync(out,Buffer.from(await bundle.blob.arrayBuffer()));
console.log(JSON.stringify({bundle:'original_MultiOmics_JS_research_ZIP_packager',synthetic_inputs_only:true,
  downstream_status:bundle.manifest.downstream_exchange.status,
  contents:bundle.manifest.downstream_exchange.matrices,
  n_zip_members:Object.keys(bundle.manifest.files).length},null,2));

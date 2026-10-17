/**
 * Self-host precisely the same variable font families as PMx Explain.
 * Runs during GitHub Actions build only; no runtime npm/CDN/tracker requests.
 * Never copy these font binaries into docs, reports, logs, attachments or chat.
 */
import { copyFileSync, existsSync, mkdirSync, statSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const site = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const destination = resolve(process.argv[2] || join(site, 'dist'));
const output = join(destination, 'assets', 'fonts');
const fonts = [
  ['inter', 'inter-latin-wght-normal.woff2'],
  ['source-serif-4', 'source-serif-4-latin-wght-normal.woff2'],
  ['jetbrains-mono', 'jetbrains-mono-latin-wght-normal.woff2'],
];
mkdirSync(output, { recursive: true });
for (const [name, filename] of fonts) {
  const source = join(site, 'node_modules', '@fontsource-variable', name, 'files', filename);
  if (!existsSync(source)) {
    throw new Error('Required original @fontsource-variable latin font is missing: ' + name);
  }
  if (statSync(source).size < 10000) {
    throw new Error('Unreasonable empty or malformed font asset: ' + name);
  }
  copyFileSync(source, join(output, filename));
}
console.log('Three pinned @fontsource-variable Latin fonts self-hosted in the static build');

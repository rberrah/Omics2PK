# Proposed license strategy (decision pending — 2026-10-16)

**Status:** advice and owner review only. **No license is added or changed by this document.** The intention under discussion is to allow free academic research and sharing, while requiring specific authorization for commercial exploitation.

## 1. Two distinct goals

- **Open source / OSI-approved**: licenses such as MIT, Apache-2.0, GPL-3.0 or AGPL-3.0 allow commercial use. AGPL requires some source-sharing obligations for modified network services, but does not prohibit paid use.
- **Source-available noncommercial**: **PolyForm Noncommercial 1.0.0** is a potential standardized software license for new, rights-cleared code. It permits use, changes and distribution for noncommercial purposes. Its explicit permitted-purpose clause includes educational institutions, public research organizations and charitable bodies regardless of funding source. It **is not an OSI open-source license**. Commercial companies may need a separately negotiated license even for research. See official text: https://polyformproject.org/licenses/noncommercial/1.0.0

Creative Commons advises against using its CC BY-NC licenses for source code: https://creativecommons.org/faq/ . CC BY-NC may be appropriate for independently owned literary/artistic content after checking compatibility, but not as the default software license.

## 2. Existing grants cannot be undone

The following accessible rberrah repositories had a root MIT LICENSE when inspected on 2026-10-16:

- rberrah/MultiOmics
- rberrah/rberrah.github.io
- rberrah/internat
- rberrah/stats

Their public code versions already granted MIT permissions, including commercial use. Replacing a LICENSE file in a later version **does not revoke those existing grants**, nor turn upstream MIT source into exclusively noncommercial software. New separately authored components may have different terms only if authors have the right to license them, and copyright notices and license compatibility are preserved. MultiOmics also declares CC BY-SA 4.0 for certain documentation/content; existing grants remain effective.

The other accessible repositories (including Omics2PK at this stage) had no root LICENSE in the inspected default branch, but that does **not** prove ownership, absence of licenses in subdirectories, or absence of third-party restrictions.

## 3. Omics2PK decision gate

Before distributing under PolyForm Noncommercial or a dual-license model:

1. Confirm intellectual-property ownership, host institution's policies and collaborators' contributions.
2. Distinguish **newly authored Omics2PK code**, incorporated upstream code (e.g. SWEET MIT), and scientific datasets. Never remove source notices.
3. Define noncommercial academic use including funded collaborations, sponsored university projects, training workshops and clinical institutional work; assess whether PolyForm's terms align with intended partners.
4. For companies, consider a **commercial license by agreement** from the rights holder, but do not promise license grants on code not exclusively controlled.
5. If maximum uptake and industry benchmarking matter more than reserving commercial rights, retain MIT/Apache-2.0 for the method and consider separately licensing a commercial hosted application instead.

## 4. Repository-wide policy

Audit **each repository** separately; do not apply the same change mechanically to research code, a public website, teaching text, art, private institutional code, and collaborative work. New source-available research repositories may use a noncommercial license after review. Existing MIT/CC BY-SA content should keep accurate notices. There is no automatically imposed uniform noncommercial license.

**Interim state for Omics2PK:** public GitHub visibility; no explicit software license selected. This is not a declaration of open-source or permission for unrestricted reuse.

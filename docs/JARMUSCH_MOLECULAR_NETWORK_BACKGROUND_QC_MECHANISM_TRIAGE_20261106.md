# Jarmusch — audit source original des ions, fonds analytiques et priorités mécanistiques
**6 novembre 2026 — document de décision scientifique, exploratoire, aucun gain PK validé.**

## Données et reproductibilité
- Source originale GNPS Jarmusch (2019 task `b3966445abad4658a3cdc63c8198a666`), 6 349 ions et matrices biologiques originales.
- Audit des voisinages GNPS et des rétentions : [CI #42137033994](https://github.com/rberrah/Omics2PK/actions/runs/42137033994) ; script `scripts/audit_original_jarmusch_ion_graph_context.py`.
- Contrôle de détection par matrice, même n10 original, et classes de blancs : [CI #42137390335](https://github.com/rberrah/Omics2PK/actions/runs/42137390335) ; script `scripts/audit_original_jarmusch_selected_ion_quality.py`.
- **AUC/CL clinique :** n10 midazolam, prédose J1/J9 ; gain apparent plasma J1 +8,2 %, urine J9 +3,8 %, intervalles descriptifs tous deux recouvrant zéro ; pas de validation indépendante.
- **Important :** ces audits utilisent le réseau déjà construit sur toute la cohorte/source. Les informations de réseau et de rétention servent à une **annotation après découverte**, pas à la sélection des variables d'un modèle inductif externe.

## Nouveaux faits issus des spectres originaux

| Ion GNPS | Prélèvement et fréquence de sélection | m/z | RTMean original (unité à confirmer dans la source) | Liens de réseau GNPS | Étiquette actuelle |
|---|---|---:|---:|---|---|
| `369` | Urine J9, 6/10 | 312,2185 | 2,9945 | 6 voisins ; 2 à Δm ≈18,009 ; aucun voisin annoté | **Priorité analytique**, aucun mécanisme assignable |
| `3018` | Urine J9, 4/10 | 513,3695 | 5,7967 | 3 voisins, l'un ≈0,004 Da | **Suspect bruit/fond** ; pas de structure |
| `4321` | Urine J9, 3/10 | 293,2105 | 3,6225 | 1 voisin Δm +15,995 Da | **Suspect bruit/fond** ; Δm compatible oxydation, non démonstratif |
| `6171` | Plasma J1, 4/10 | 796,6213 | 5,6598 | 0 lien au seuil GNPS choisi | **À identifier** ; isolement spectral |
| `936` | Plasma J1, 2/10 | 288,1176 | 0,3706 | 2 voisins, aucun annoté | **À identifier** ; signal de rétention précoce à vérifier |
| `740` | Urine J1, 4/10 | 357,2784 | 3,6895 | Cosinus 1 avec 1889 ; même m/z | **Redondance probable** de 1889, identité présomptive CDCA |
| `1889` | Urine J1, 4/10 | 357,2785 | 3,6893 | Cosinus 1 avec 740 | **Ne pas compter comme métabolite distinct confirmé** |

**Aucune attribution de molécule ne peut reposer seulement sur m/z, Δm, rétention et réseau de similitude.** Les ions 740/1889 doivent être rassemblés pour un contrôle de redondance (relecture chromatogrammes, standards, extraction et spectres MS/MS), pas fusionnés sans vérification du brut.

## Le contrôle qualité modifie la hiérarchie

La table source permet de compter les détections dans **153 blancs d'analyse** et **5 blancs d'extraction**, dont la comparabilité analytique par lot n'est pas encore validée :

| Ion | Dans les 10 sujets, prélèvement prioritaire | Blancs d'analyse positifs | Blancs d'extraction positifs | Décision prudente |
|---|---:|---:|---:|---|
| `369`, urine J9 | 10/10 | 11/153 | 5/5 | Pas immédiatement éliminé ; confirmer le signal au-dessus du fond et dans le même lot |
| `3018`, urine J9 | 6/10 | **120/153** | 4/5 | **Alerte fond élevé**, exclure d'un panel confirmatoire tant que non qualifié |
| `4321`, urine J9 | 10/10 | **131/153** | 5/5 | **Alerte fond élevé**, ne pas attribuer à CYP3A |
| `6171`, plasma J1 | 9/10 | 54/153 | 3/5 | Alerte modérée ; contrôler blanc/échantillon et lot |
| `936`, plasma J1 | 7/10 | 6/153 | 5/5 | Faible détection dans les blancs d'analyse, mais blancs d'extraction positifs |
| `6135`, urine J9 | 10/10 | 136/153 | 5/5 | Red flag fond/dérive; présence dans composante biliaire GNPS insuffisante |

**Ne jamais conclure « contaminant » sur fréquence seule.** Prochaine décision fondée sur la distribution des *intensités* d'échantillons et blancs **comparables par lot/date** et une estimation du rapport signal/fond; détecter dans un blanc n'équivaut pas à une quantité significative. Les prélèvements prédose plasma du jour1 ont été analysés sur 5 dates distinctes et ceux du jour9 sur 6 dates; l'urine J1 et J9 sur 2 et 3 dates respectivement : risque de confusion temporelle/lot à documenter.

Ce contrôle empêche de promouvoir le petit gain PK du panel 6 349 ions comme mécanisme CYP3A sans vérifier d'abord sa composante analytique. Toute sélection de nouveaux filtres sur les mêmes dix sujets est **post-hoc exploratoire** ; elle ne constitue pas de validation externe.

## Nouveau prior mécanistique de littérature (2026)

**Yee et al. 23 septembre 2026, CPT, DOI [10.1002/cpt.70489](https://doi.org/10.1002/cpt.70489)** : chez 11 volontaires exposés à l'itraconazole, la mesure ciblée de GHCA, GUDCA et conjugués de 1β-hydroxy-deoxycholate présente des diminutions marquées sous inhibition CYP3A. Les auteurs caractérisent *in vitro* l'hydroxylation **CYP3A4-sélective** du glycochenodeoxycholate (GCDCA) vers la glycohyocholate (GHCA). C'est une hypothèse biochimique prioritaire, **pas une preuve que GHCA prédit la PK de tacrolimus, ni un marqueur CYP3A5 spécifique**. L'étude cible une inhibition CYP3A4 avec un inhibiteur, pas le même mécanisme que l'antibiothérapie de Jarmusch; elle ne réplique pas les n10.

**Rodrigues 2026, J Clin Pharmacol, DOI [10.1002/jcph.70165](https://doi.org/10.1002/jcph.70165)** : orale midazolam combine composantes hépatique et intestinale ; les marqueurs plasma/urine peuvent refléter différemment CYP3A hépatique et intestinal. Variations genotype CYP3A5, transporteurs, rythme circadien, élimination rénale et assay time doivent être explicitement modélisés.

**Biomarqueurs à tester biologiquement AVANT estimation PK** :
1. Plasma ciblé GCDCA et **GHCA** (substrat/produit CYP3A4, avec identification et valeurs absolues) ; GUDCA comme marqueur potentiellement moins spécifique.
2. Plasma/urine 1β-hydroxy-deoxycholate / deoxycholate et conjugués glycine/taurine avec contrôle d'excrétion et microbiote.
3. Urine cortisol/6β-hydroxycortisol, cortisone/6β-hydroxycortisone, **rapport mesuré** et prélèvement circadien standardisé.
4. Plasma 4β-hydroxycholestérol/cholestérol, cinétique lente et possible faible pouvoir prédictif à court terme.
5. Pour CYP3A5*3, la preuve proximale reste l'**épissage aberrant** et la baisse de protéine active dans le foie/intestin ; les marqueurs plasmatiques CYP3A4 ne peuvent pas remplacer un test spécifique CYP3A5.

**Trois règles de qualification** : (A) bon organe / bon prélèvement / bon temps ; (B) identité et réaction biochimique vérifiées ; (C) gain incrémental PK évalué indépendamment de la sélection et des données de découverte. La décision de k=0 est autorisée à l'étape C, mais ne dispense pas de vérifier A et B.

**Décision actuelle :** poursuivre en priorité une qualification spectrale et de fond pour 369/936/6171, ne pas retenir 3018/4321 comme acteurs biologiques sans requalification, et conserver les hypothèses GHCA/GCDCA et 1βOH-DCA/DCA pour un véritable dosage ciblé prospective/source indépendant.
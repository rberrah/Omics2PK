# Omics2PK — Lawson / Jardou IATDMCT 2025 : protocole de validation mécanistique et PK

**État au 6 novembre 2026 : données en cours d'obtention auprès de l'équipe ; aucun fichier individuel reçu ou exploité ici. Ce document ne prétend pas reproduire les résultats.**

## Source originale (identifiée)

- IATDMCT 2025, Singapour, résumé **326**, « Lipopolysaccharide-linked endotoxemia correlates with cyclosporine A exposure in transplant patients » ; Manon Jardou, Clarisse Brossier, Alexandre Garnier, Chloé Barny, Hélène Roussel, Antoine Humeau, Emilie Pinault, François-Ludovic Sauvage, Caroline Monchaud, Pierre Marquet, Anne Druilhe, Roland Lawson.
- Source primaire : https://www.iatdmct2025.org/wp-content/uploads/2025/09/IATDMCT2025-Abstract-book.pdf
- 48 transplantés : PIGREC **cœur**, STIMMUGREP **poumon**, BIOMARGIN **rein** ; prélèvements conservés à CRBIOLIM à **un an post-greffe**. Témoins : 24 volontaires sains, **sans AUC immunosuppresseur comparatrice**.
- **Biomarqueur :** lipopolysaccharide **LPS plasmatique**, dosage **ELISA**, constituant d'endotoxémie lié aux bactéries Gram négatif et à la perméabilité intestinale. **C'est une mesure ciblée mono-analyte, pas un transcriptome, protéome, métabolome, métagénome ou une mesure directe de β-glucuronidase.**
- **PK :** AUC individuelles de ciclosporine (CsA), tacrolimus et mycophénolate (MMF administré, **analyte PK à vérifier : MPA ?**), extraites de bases historiques. Résultat du résumé : seule **ciclosporine AUC** corrélée au LPS (ρ de Spearman **0,818**, p=**0,002**). **Les effectifs CsA/TAC/MPA par bras ne sont pas rapportés** ; ne jamais présenter ρ comme calculé sur les 48 patients.
- Le résumé indique les moyennes LPS 99±22 pg/mL et 257±32 pg/mL (SEM), et dit explicitement LPS plus élevé chez les transplantés. **L'ordre de présentation des groupes dans la phrase source est ambigu** ; vérifier tableau primaire avant réutilisation chiffrée.

## Ce que la cohorte peut tester, et ce qu'elle NE peut pas tester

- **E3 clinique, incrément prédictif à démontrer** : ajouter **LPS prédose** à un modèle PK clinique/dose et CYP/transporteur PGx si les génotypes et temps existent. Une **corrélation marginale ρ** n'est ni un gain de prédiction, ni un effet causal.
- La greffe et le médicament sont **colinéaires partiellement** : cœur/poumon/rein peuvent avoir des protocoles CsA/TAC différents. Ne pas ajuster naïvement un effet « CsA » en poolant trois organes avec peu de chevauchement.
- **Temporalité :** si LPS a été mesuré en même temps ou après l'AUC, le résultat est un **biomarqueur concomitant**, pas un prédicteur prospective d'exposition. Exiger date/heure/aliquote et fenêtre d'AUC.
- **Causalité réciproque :** LPS pourrait moduler inflammation/CYP/transporteurs → PK ; les immunosuppresseurs et MMF pourraient modifier le microbiote, la perméabilité intestinale et le LPS → rétroaction. Risque de confusion (infection, inflammation, anti-infectieux, corticothérapie, dysbiose, fonction hépatique).
- **MPA mécanistique distincte de LPS :** **MPAG (biliaire) → activité bactérienne β-glucuronidase intestinale → MPA réactivé → recirculation entérohépatique et toxicité digestive**. Un LPS plasmatique seul ne mesure **ni β-glucuronidase, ni microbiome, ni le pic secondaire**. Source Lawson/Brossier https://doi.org/10.1016/j.lfs.2024.122792.
- **Toute greffe rénale :** distinguer le CYP3A5 du receveur pertinent pour la PK systémique (foie/intestin) et celui du donneur pertinent surtout localement. Pour **greffe hépatique**, il faudrait distinguer donneur foie et receveur intestin (non inclus dans les 3 sous-cohortes présentées).
- Les témoins sains peuvent vérifier l'endotoxémie/qualité de biomarqueur mais **ne doivent pas devenir des observations pour l'entraînement d'AUC inexistantes**.

## Contrat de données attendu : ne rien inventer avant réception

| Champ nécessaire | Pourquoi il est indispensable | Niveau actuel |
|---|---|---|
| identifiant pseudonymisé **même personne** pour LPS, AUC, éventuelles omiques/PGx | lien patient prouvé, prévenir fausses jointures | NON reçu |
| origine PIGREC/STIMMUGREP/BIOMARGIN, organe greffé, **type de médicament** | stratifier co-localisation des protocoles | résumé seul |
| dose, formulation, date/heure dernière dose, profil de prélèvements, **méthode de calcul de l'AUC et fenêtre** | éviter AUC non comparable, contrôle exposition et erreurs LSS | NON reçu |
| LPS ELISA pg/mL, date et heure de prélèvement, répétitions techniques, stockage, hémolyse | bioanalytique, antériorité, effet de lot | résumé seul |
| génotype clinique réel CYP3A4*22, CYP3A5*3/*6/*7, ABCB1, UGTs **si disponible** | comparateur PGx complet ou partiel ; ne jamais inférer variant absent | disponibilité NON confirmée |
| Ht, CRP, enzymes hépatiques, albumine, créatinine/eGFR, antibiotiques/azolés, stéroïdes, temps postgreffe | facteurs de confusion et mécanisme contextuel | NON reçu |
| données MPA/MPAG *simultanées* si disponibles ; profil après dose | relier éventuelle activité β-glucuronidase à la composante entero-hépatique | NON confirmé |
| profils de microbiote/fèces, β-glucuronidase, métabolites identifiés, transcriptome/protéome | vraie pharmacomicrobiomique multimodale plutôt que LPS seul | NON confirmé |
| consentements/autorisations de réutilisation et partage de pseudonymisation | reproductibilité sans diffusion des dossiers patients | à vérifier via gouvernance |

### Gate d'admissibilité dans Omics2PK

- `E3-LPS` si LPS plasmique + vraie AUC du **même sujet** ; `E3-omics` uniquement si d'autres mesures moléculaires biologiques identifiées avec tissu/temps sont appariées ; `E3-PGx-increment` uniquement si génotypes individuels réels disponibles.
- AUC **CsA** d'abord, car seule association publiée ; TAC et MPA comme contrôles biologiques et médicaments différents avec **hypothèses distinctes**, pas comme réplicats d'un même effet.
- La date LPS vs AUC doit être décrite **avant** de qualifier le design comme « prédictif » ; sans antériorité prouvée, ne rapporter que l'association/prediction concomitante exploratoire.
- Taille par médicament **n_CsA / n_TAC / n_MPA** et nombre par organe doivent précéder toute CV. Le n=48 global ne prouve pas que chaque bras est statistiquement estimable.

## Plan d'analyse préspécifié (à verrouiller une fois les données auditées)

**A. Plausibilité mécanistique** : vérifier un lien LPS ↔ phénotype inflammatoire / hépatique et microbiome, si ces mesures sont réellement présentes ; pour MPA, regarder ratio MPA/MPAG et pic secondaire, β-glucuronidase si mesurée. Chercher cohérence de direction **sans imposer** l'effet.

**B. Association clinique** : reproduire Spearman ρ, p et CI de CsA, TAC, MPA sur les effectifs réels ; tracés LPS vs dose-normalized AUC avec erreurs d'AUC et ajustement robuste selon taille (dose, organe, CRP et inhibiteurs si disponibles). Rapporter p descriptif et correction pour au moins 3 médicaments.

**C. Valeur prédictive** : pour **chaque médicament séparément**, modèles incrémentaux `clinical` → `clinical+PGx` → `clinical+PGx+LPS`, éventuel `+B-glucuronidase / omiques documentées`. Baseline contient dose et covariables antérieures, le comparateur ne doit pas inclure une **concentration médicamenteuse issue du profil utilisé pour calculer l'AUC**, sinon fuite par construction. LOPO avec validation interne pour pénalisation, groupes par sujet ; bloc de régularisation/parcimonie et k0 autorisé. MAE log AUC, erreur relative, P20, biais et CI appariées descriptives. En cas de trop faible n par bras, **NO-GO model selection** : association exploratoire uniquement, ou modèle à 1 covariable a priori sans tuning multiple.

**D. Transfert causal** : comparer direction des effets selon organe et co-médicaments ; ne jamais déclarer causal à partir d'un ρ transversal. En cas de données répétées, stratégies **within-patient** permettant de réduire la confusion interpersonnelle.

**E. Article/communication** : *plasma endotoxemia as a clinically accessible but mechanistically indirect biomarker; prespecified mechanistic test and incremental PK utility over known covariates*. Faire apparaître la piste β-glucuronidase/MPA comme une **hypothèse à tester** plutôt qu'un résultat Lawson déjà démontré.

## Études publiques en complément (et pourquoi)

- **Sarömba et al. 2026 n47** : plasma prédose prochaine dose, CYP2D6 génotypes et deux rapports métaboliques nommés ; vraie concentration métoprolol sous traitement ; Omics2PK reproduit 25,94% de baisse MAE avec les deux rapports versus clinique+PGx sur original n47 en LOPO (pas AUC, sans cohorte externe). https://github.com/rberrah/Omics2PK/actions/runs/42179780471
- **Müller et al. 2026 n20** : même chimie CYP2D6, plasma après repas pomme de terre standardisé à 9 temps 0–29 h ; vérifie alimentaire/temps sur le phénotypage plutôt que gain PK ; original XLSX disponible https://github.com/JulianPeterMueller/Solanidine-standardized-potato-meal.
- **Rodrigues et al. 2021** : protéomique ciblée **CYP3A4/5 des vésicules extracellulaires d'origine hépatique dans sérum**, vrais sujets avec induction rifampicine et AUC midazolam oral correspondante (groupes n4 et n5), génotype CYP3A5*1/*3 chez les sujets signalés ; groupe très petit, valeurs patient-à-patient non qualifiées. https://doi.org/10.1002/cpt.2244
- **Newman et al. 2024** : vrai foie humain n11, **homogénat ↔ vésicules extracellulaires isolées du foie**, même donneur, protéines CYP3A4 et CYP3A5 ciblées ; corrélation CYP3A4 ρ=0,791. **Attention, EV du tissu hépatique ne signifie pas EV du plasma sanguin.** https://doi.org/10.1111/bcp.16187
- **Lawson 2025** : transferts rein/cœur/poumon + LPS plasmique + véritables AUC, mais **génétique non établie** dans le résumé; unique cohorte de transplantation de ce groupe, ne pas confondre avec jeux de données publics.

**Décision** : ne pas attribuer de nouveau gain Omics2PK aux 48 Lawson tant que source individuelle et disponibilité PGx/omiques ne sont pas vérifiées. Le protocole et ses gates sont prêts pour les agents.

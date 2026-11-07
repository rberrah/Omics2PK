# ARIANA DRKS00025716 — qualification de la forme galénique du métoprolol (7 novembre 2026)

**Conclusion scientifique : données fortement compatibles avec du métoprolol succinate à libération prolongée ; la formulation exacte au niveau patient n'est PAS confirmée.**

## Trois niveaux de preuve (ne pas fusionner)

1. **Enregistrement d'essai source (vérifié)** : German Clinical Trials Register [DRKS00025716](https://drks.de/search/de/trial/DRKS00025716). ARIANA = `All-comer` Registry for ImmunocArdiology aNd cardiometabolic disease Aachen. Inscription 30 août 2021, début réel 7 février 2022, *registre prospectif non interventionnel*. Il ne prescrit **pas** une spécialité de métoprolol et indique **aucun protocole ou document complémentaire déposé** (rubrique `Studienprotokoll und weitere Studiendokumente: Kein Eintrag`). On ne peut donc pas déduire la formulation du registre.
2. **Article publié — primaire (vérifié)** : Sarömba et al., *Frontiers in Pharmacology* 2026, [DOI 10.3389/fphar.2026.1812007](https://doi.org/10.3389/fphar.2026.1812007), Methods/Results. Les 47 patients analysés avaient un traitement chronique stable ≥1 semaine, **deux prises égales le matin et le soir**, et une dose **journalière** de **95 mg (n=23)** ou **190 mg (n=24)**. Prélèvement à jeun avant la prise suivante (Cmin supposée) ; ni marque commerciale, ni tartrate/succinate, ni mention explicite `retard`/IR/OAD dans l'article. Les heures exactes dose/prélèvement ne sont pas communiquées.
3. **Notice réglementaire de spécialités allemandes (correspondance indirecte forte)** : [Metoprolol succinate TAD — Rote Liste / notice allemande](https://www.patienteninfo-service.de/a-z-liste/m/metoprolol-succinat-tadR-2375-mg-475-mg-95-mg-190-mg-retardtabletten), ou [Metoprololsuccinat Heumann 2021](https://www.heumann.de/fileadmin/user_upload/produkte/infos/Fachinformation-Metoprololsuccinat-Heumann-47_5-95-190-mg-Retardtabletten.pdf). Dosages des **comprimés LP de succinate** : 47,5 mg, 95 mg, 190 mg ; ces doses sont équivalentes respectivement à **50 mg, 100 mg, 200 mg de tartrate** en teneur active. Schémas ARIANA :
   - 95 mg/j en deux prises -> **47,5 mg matin + 47,5 mg soir**, exactement un dosage LP commercialisé ;
   - 190 mg/j en deux prises -> **95 mg matin + 95 mg soir**, exactement un autre dosage LP.
   - Les notices LP recommandent **habituellement une administration quotidienne** ; ici les auteurs décrivent deux administrations quotidiennes. Le BID **n'exclut pas la LP**, mais interdit de supposer que le schéma suit l'AMM standard du fabricant.

## Classification

- `formulation_confirmed`: **false**.
- `published_dose_strength_consistent_with_succinate_retard`: **true**.
- `formulation_working_hypothesis`: **succinate_LP_BID**, avec statut **fortement suggérée, pas vérifiée**.
- `alternative`: tartrate immédiat et/ou dossier clinique posologique exprimé en équivalent succinate, **possible mais moins étayé**.
- **Attention unité :** 95 mg de succinate LP = équivalent nominal **100 mg de tartrate**. Confirmer convention de dose du modèle PopPK (masse de sel, de métoprolol base ou équivalent tartrate) avant simulation ; éviter l'erreur systématique de dose.
- **Horaires manquants :** ne pas remplacer `prédose` par `t=12 h exact` sans analyse de sensibilité. La LP BID et l'IR BID produisent des concentrations pré-dose différentes, même avec doses nominales équivalentes.

## Conséquences pour l'évaluation Omics2PK / ARIANA v3

1. Mettre **LP succinate BID** en scénario pharmaceutique principal *hypothétique* et IR BID en scénario de sensibilité ; préciser la source et maintenir un drapeau `formulation_unverified`.
2. Pour le comparateur PopPK, faire valoir un modèle **LP avec fonction d'absorption adaptée** (de préférence modèle publié avec vraie LP) ; Eugene/Chai IR sont des contrôles secondaires et ne doivent pas être présentés comme modèles structurels confirmés d'ARIANA.
3. Comparer des concentrations pré-dose aux modèles sur une grille d'intervalles plausibles et dose préalable documentée ; aucune estimation MAPBE utilisant le Cmin cible.
4. Ne jamais conclure à la validité d'une galénique en choisissant celle qui donne la meilleure MAE sur les 47 personnes et en réévaluant ces mêmes sujets : ce serait un **choix post-hoc** optimisé sur le test.
5. Demander aux auteurs : **sel utilisé, marque commerciale, nom exact de spécialité + code ATC/PZN si disponible, IR/LP, fractionnement BID, durée de traitement, heures des deux doses et du prélèvement, unités dose** ; vérifier si homogène sur les 47 patients.
6. Contacts : correspondant de l'article **Katja S. Just** (kjust@ukaachen.de) ; contact scientifique de l'enregistrement **Berkan Kurt** (bkurt@ukaachen.de). Aucun contact n'a été envoyé.

**Source de vérification sur registre** : https://drks.de/search/de/trial/DRKS00025716 (rubriques `Studiendesign`, `Studienprotokoll`, `Weitere Studiendokumente`). **Source de l'article** : https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2026.1812007/full . **Source posologies LP** : https://www.patienteninfo-service.de/a-z-liste/m/metoprolol-succinat-tadR-2375-mg-475-mg-95-mg-190-mg-retardtabletten .

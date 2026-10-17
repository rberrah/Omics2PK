# 17 octobre 2026 — recherche de cohortes cliniques indépendantes : metformine (n17) et contrôles négatifs

## Nouvelle piste prioritaire : Tee et al. 2022, metformine PK–pharmacometabolomique

**Article principal authentifié :** Tee et al., *Pharmaceutics* 2022, DOI https://doi.org/10.3390/pharmaceutics14061268 , PMID 35745841, [PMC9231303](https://pmc.ncbi.nlm.nih.gov/articles/PMC9231303/), **CC BY** et supplément annoncé chez l'éditeur https://www.mdpi.com/article/10.3390/pharmaceutics14061268/s1 . Source audit dans le dépôt : `scripts/audit_metformin_17_pharmacometabolomic_public.py` ; [CI original 17 octobre](https://github.com/rberrah/Omics2PK/actions/runs/39019261323).

- **17 volontaires sains**, dose orale réelle de **1 000 mg de metformine**, **15 temps de prélèvement plasma PK par personne**, soit **255 mesures de concentration** présentées comme quantifiées. Article et légendes identifient bien les valeurs PK spécifiques aux volontaires ; supplément **Table S5.1** : concentrations selon les temps, **S5.2** : paramètres PK individuels.
- Expérience **métabolomique LC-MS non ciblée du plasma avant la dose** et aux temps guidés par le Tmax/Cmax individuel ; certains ensembles contiennent 17 paires pré-dose/maximum et d'autres seulement 6 à 18 paires selon le sous-protocole. Ce protocole offre potentiellement un vrai exemple *pré-dose metabolome→future PK*, plus utile pour E3 que la plupart des études transversales de greffe.
- **Point critique :** la publication montre surtout PCA/PLS-DA, pics/features filtrés et voies moléculaires enrichies ; aucune **matrice brute par sujet de features métabolomiques appariée** à la table des 255 concentrations n'a été vérifiée comme téléchargeable. Les seules signatures de voies groupées ou listes de métabolites ne remplacent pas la matrice individuelle.
- L'audit réel a récupéré **les métadonnées officielles PMC OA** du manuscrit `PMC9231303.1`, licence **CC BY**, mais **zéro objet supplémentaire Office** dans la liste des fichiers stockés en libre accès au sein de PMC OA. Cela ne prouve pas que le supplément de l'éditeur soit absent : il est référencé par la publication mais **le script ne l'a pas extrait de l'API officielle testée**. Le score `paired_PK_metabolomics_available=False` est volontaire.

### Prochaine qualification binaire (E3)

1. Vérifier, dans **Table S5 originale de l'éditeur**, que les 17 sujets ont des IDs et leurs **15 temps, unités ng/mL, traitement BLQ, dose 1 000 mg, temps relatif** ; le lien « supplément disponible » ne démontre pas encore la qualité du fichier binaire.
2. Vérifier la disponibilité d'une **vraie matrice de m/z/intensités métabolomiques pré-dose à n=17 personnes**, leur ordre d'identité, batch, corrections et échantillons de contrôle. Si seulement 6 paires pour dataset A ou matrices de PCA, interdire toute reconstruction par histogrammes.
3. Vérifier `subject_id` **strictement identique** PK↔MS et que la signature moléculaire utilisée pour prédire la PK a été prélevée **avant la dose**, non après Cmax (sinon fuite inverse / marqueur de réponse au médicament).
4. En cas de données individuelles originales ouvertes et légalement utilisables, évaluer sur donneur-exclu `baseline clinical (sex, age, weight)` vs `pre-dose LC-MS sparse regularized scores`, sur les **17 mêmes individus**, sortie Cmax/AUC et NCA vérifiée à partir des 255 concentrations, sans chercher le modèle de grande dimension optimal sur ces 17.
5. Ne jamais déclarer E3 validé tant que l'étape **métabolites pré-dose × patient × C(t)** ne passe pas.

## Contrôle indépendant tacrolimus/16S : Degraeve n93

- Article DOI https://doi.org/10.1002/cpt.3077 et SRA BioProject https://www.ncbi.nlm.nih.gov/bioproject/PRJNA929904 (93 entrées de séquençage). La cohorte comprend **93 transplantés rénaux CYP3A5 non-exprimeurs**, original TAC C0 et doses cliniques plus estimations d'AUC basées sur limited sampling : c'est une source cliniquement très proche de la thèse.
- Le supplément original **avait déjà été audité** et ne contient aucune table vérifiée de 93 identifiants patient↔SRA/ASV↔C0/dose↔C1/C3/AUC. Les objets GraphPad et tableaux par groupes **ne constituent pas** des matrices individuelles prêtes à prédire. [Audit antérieur du vrai ZIP](DEGRAEVE_WILEY_UPLOADED_SUPPLEMENT_AUDIT_20261016.md).
- La présence d'exactement 93 BioSamples dans SRA ne prouve **pas** la possibilité de relier PK par identifiant au niveau de chaque individu. Il ne faut jamais joindre « par ordre de ligne ».

## Autres études : état et priorité

| Source | Molécule | Omique et PK observées dans même protocole | PK originelle numérique libre | Matrice moléculaire↔PK individuelle vérifiée | Décision |
|---|---|---|---|---|---|
| Tee 2022, n17 | Metformine | Plasma metabolomics baseline + dose 1000 mg + C(t) 15 temps | Annoncée en supplément Table S5.1, extraction à qualifier | **Non** | **Nouveau P0 : priorité à la donnée individuelle des deux matrices** |
| Degraeve 2024, n93 | Tacrolimus | 16S fécal + C0/dose + AUC estimée | Seulement groupe/régression confirmé, original sujet↔SRA non | **Non** | P1 : acquisition source clinique + ID croisé, pas de modèle forcé |
| Zhu 2022, n137 / 693 C0 | Tacrolimus | Métabolomique sang + TDM longitudinal | Non vérifiée publiquement | **Non** | P1 : meilleure cible E4 mais probablement demande de données |
| Tee 2021, n26 | Metformine | Métabolomique longitudinale et Cmax individuels mesurés | Concentrations PK dites **unpublished data** | **Non** | Ne pas interpréter le suivi métabolomique comme PK disponible |
| Phapale 2010, n29 | Tacrolimus | Métabolomique urine pré-dose + PK 0–72 h | Profil PK présenté dans publication ; 13 suppléments audités sans appariement source original complet | **Non** | P1, original complet à demander |
| NCI60, 60 lignées | Doxorubicine | Multiomics et GI50 | GI50 est réponse et cytotoxicité, pas PK | Non applicable | Exclure d'une validation Omics→PK patient |

**Bilan scientifique :** les liens biologiques omiques→activité CYP/NAT ont des données originales vérifiées ; **aucun jeu clinique totalement ouvert et joint Omics pré-dose × mesures plasmatiques dose/temps/C(t) n'est encore qualifié**. Cette affirmation conserve sa force jusqu'à vérification de tous les fichiers individuels, pas seulement de l'existence d'une publication.


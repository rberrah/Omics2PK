"""Render the Omics2PK static research portal, same editorial design as PMx Explain.

No JavaScript framework, tracking, external fonts, data upload, or clinical claims.
Static project Pages paths: /Omics2PK/{demo,guide,analyze,ecosystem,disciplines,resources}/
The legacy /learn/ route is retained as a compatibility alias for the guided course.
"""
from __future__ import annotations
import argparse
from html import escape
from pathlib import Path

SITE_TITLE = "Omics2PK"
PMX = "https://rberrah.github.io/pharmacometrie/"
MAIN = "https://rberrah.github.io/"
MULTIOMICS = "https://rberrah.github.io/multiomics/"
REPO = "https://github.com/rberrah/Omics2PK"
ROOT = "https://rberrah.github.io/Omics2PK/"

def t(fr,en, tag="span", classes=""):
    attr = (" class='"+escape(classes,quote=True)+"'") if classes else ""
    return f'<{tag}{attr} data-fr="{escape(fr,quote=True)}" data-en="{escape(en,quote=True)}">{escape(fr)}</{tag}>'

def a(path,label_fr,label_en,*,base,external=False):
    href = path if external else base+path
    ext=' target="_blank" rel="noopener noreferrer"' if external else ""
    return f'<a href="{escape(href,quote=True)}"{ext} data-fr="{escape(label_fr,quote=True)}" data-en="{escape(label_en,quote=True)}">{escape(label_fr)}</a>'


# Inline original Lucide icon SVG paths, same glyphs and sizes used by the Svelte sites.
# Decorative icons contain no interactive handler and do not require a client bundle.
LUCIDE = {
    "demo": '<path d="M9 3h6"/><path d="M10 9h4"/><path d="M10 3v6l-5 9a2 2 0 0 0 1.7 3h10.6a2 2 0 0 0 1.7-3l-5-9V3"/><path d="M7.5 15h9"/>',
    "guide": '<path d="M12 5v16"/><path d="M20.001 19A2 2 0 0022 17V5a2 2 0 00-1.999-2L16 3.002A5 5 0 0012 5a5 5 0 00-4-2H4a2 2 0 00-2 2v12a2 2 0 001.999 2H8a5 5 0 014 2 5 5 0 014-2z"/>',
    "analyze": '<path d="M12 16v5"/><path d="M16 14.639V21"/><path d="M20 10.656V21"/><path d="m22 3-8.646 8.646a.5.5 0 0 1-.708 0L9.354 8.354a.5.5 0 0 0-.707 0L2 15"/><path d="M4 18.463V21"/><path d="M8 14.656V21"/>',
    "disciplines": '<path d="m16 6 4 14"/><path d="M12 6v14"/><path d="M8 8v12"/><path d="M4 4v16"/>',
    "resources": '<path d="m16 6 4 14"/><path d="M12 6v14"/><path d="M8 8v12"/><path d="M4 4v16"/>',
    "down": '<path d="m6 9 6 6 6-6"/>',
}

def lucide(name, size=16):
    return (f'<svg class="lucide-icon" width="{size}" height="{size}" '
            'viewBox="0 0 24 24" fill="none" stroke="currentColor" '
            'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" '
            'aria-hidden="true">' + LUCIDE[name] + '</svg>')

def nav(base,page):
    groups = [
      ("demo","Pourquoi ?","Why?",[
         ("demo/","Démo scientifique","Scientific demo",False),
         ("demo/#conclusion","Ce que les preuves disent","What the evidence says",False)
      ]),
      ("guide","Comment faire ?","How?",[
         ("guide/","Cours guidé","Guided course",False),
         ("guide/#checklist","Checklist d’une étude","Study checklist",False)
      ]),
      ("analyze","Analyser","Analyze",[
         ("analyze/","En cours de développement","In development",False)
      ]),
      ("disciplines","Écosystème","Ecosystem",[
         ("ecosystem/","Vue d’ensemble","Ecosystem overview",False),
         (PMX,"PMx Explain ↗","PMx Explain ↗",True),
         (MULTIOMICS,"MultiOmics ↗","MultiOmics ↗",True)
      ]),
      ("resources","Ressources","Resources",[
         ("resources/","Références et documentation","References and documentation",False),
         (REPO,"Code source et rapports ↗","Source code and reports ↗",True),
         (MAIN,"Site principal ↗","Main website ↗",True)
      ])
    ]
    details=[]
    for key,fr,en,links in groups:
        inner="".join(a(path,lfr,len_,base=base,external=ext) for path,lfr,len_,ext in links)
        details.append('<details class="goal" name="site-objective" data-testid="goal-'+key+'">'
          '<summary'+(' class="active"' if (page==key or (key=="guide" and page=="learn") or (key=="disciplines" and page=="ecosystem")) else '')+'>'+lucide(key,16)
          +t(fr,en)+lucide("down",13)+'</summary><div class="goal-links">'+inner+'</div></details>')
    return ('<header class="site-header"><a class="logo" href="'+base+'" aria-label="Omics2PK, accueil" data-testid="logo-link">'
        '<span class="mark" aria-hidden="true">oP</span><span class="word">Omics<em>2PK</em></span></a>'
        '<button class="burger" type="button" aria-label="Ouvrir le menu" aria-controls="site-nav" aria-expanded="false" data-testid="nav-toggle"><span></span><span></span><span></span></button>'
        '<nav id="site-nav" aria-label="Navigation principale" data-testid="primary-nav">'
        +''.join(details)
        +'<div class="tools"><div class="language-toggle" aria-label="Language" data-testid="language-toggle">'
        '<button type="button" data-set-lang="en" aria-pressed="false">EN</button>'
        '<button type="button" data-set-lang="fr" aria-pressed="true">FR</button></div>'
        '<button type="button" class="theme-toggle" id="theme-switch" title="Thème : système" aria-label="Thème : système" data-testid="theme-toggle">'
        '<svg viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round">'
        '<rect x="3" y="4" width="18" height="13" rx="2"/><path d="M8 20h8"/></svg></button></div></nav></header>')

def card(kicker,title_fr,title_en,desc_fr,desc_en,href,link_fr,link_en,*,base):
    return ('<article class="card"><div class="card-kicker">'+escape(kicker)+'</div>'
        +t(title_fr,title_en,'h3')+t(desc_fr,desc_en,'p')
        +'<div class="card-link">'+a(href,link_fr+' →',link_en+' →',base=base,external=href.startswith('http'))+'</div></article>')

def home(base):
    return (
        '<section class="hero"><p class="eyebrow">Omics2PK · research laboratory</p>'
        +t("Des omiques à la PK : une hypothèse à tester, pas une promesse.",
           "From omics to PK: a hypothesis to test, not a promise.",'h1')
        +t("Omics2PK demande si un état moléculaire mesuré chez un individu apporte une information pharmacocinétique au-delà des covariables cliniques et pharmacogénétiques. Le site sépare volontairement la question scientifique, la méthode pour y répondre et le futur outil d’analyse.",
           "Omics2PK asks whether an individual's measured molecular state adds pharmacokinetic information beyond clinical and pharmacogenetic covariates. The site deliberately separates the scientific question, the method used to test it, and the future analysis tool.",'p','lede')
        +'<div class="cta">'+a("demo/","Voir pourquoi","See why",base=base)
        +'<a class="btn btn-outline" href="'+base+'guide/" data-fr="Apprendre comment faire" data-en="Learn how to do it">Apprendre comment faire</a></div>'
        +t("Recherche exploratoire — valeur ajoutée des omiques en PK non démontrée de façon générale ; aucun modèle de dose clinique validé.",
           "Exploratory research — general added value of omics for PK is not established; no clinically validated dosing model.",'p','hero-note')+'</section>'
        +'<section class="objectives" aria-labelledby="objective-title">'+t("Par où commencer ?", "Where should I start?",'h2')+
        '<div class="objectives-grid">'+
        a('demo/','Pourquoi tester les omiques en PK ?','Why test omics in PK?',base=base)+
        a('guide/','Comment concevoir l’étude ?','How should the study be designed?',base=base)+
        a('analyze/','Analyze — en cours','Analyze — in development',base=base)+
        a('ecosystem/','Voir l’écosystème','Explore the ecosystem',base=base)+
        '</div></section><section class="section"><p class="eyebrow">Scientific status</p>'
        +t('Ce que les données actuelles permettent déjà de dire', 'What current evidence already allows us to say','h2')
        +'<div class="cards">'
        +card('WHY · biological plausibility','ARN hépatique et activité CYP','Liver RNA and CYP activity','Sur 100 foies humains, certains transcrits CYP ajoutent un signal prédictif de l’activité enzymatique au-delà de l’âge et du sexe. Cela justifie la question biologique, pas une conclusion de PK clinique.','In 100 human livers, some CYP transcripts add predictive signal for enzyme activity beyond age and sex. This supports the biological question, not a clinical PK conclusion.','demo/#epa','Comprendre pourquoi','Understand why',base=base)
        +card('WHY · phenotype vs genotype','État fonctionnel et génotype','Functional state and genotype','Chez 13 volontaires CYP3A, un phénotype fonctionnel antérieur par sonde midazolam prédit mieux la PK ultérieure que le génotype dans ce petit jeu. Ce phénotype n’est toutefois pas une omique.','In 13 CYP3A volunteers, a prior midazolam functional phenotype predicts later PK better than genotype in this small dataset. That phenotype is not an omics measurement.','demo/#jarmusch','Voir la nuance','See the nuance',base=base)
        +card('WHY · clinical reality','Métabolome pré-dose → PK metformine','Predose metabolome → metformin PK','Premier benchmark humain Omics→PK authentique du dépôt : 32 personnes avec métabolomique pré-dose et AUC/CLr mesurées. L’ajout des omiques au modèle clinique n’améliore pas clairement la prédiction.','First authentic human Omics→PK benchmark in the repository: 32 people with predose metabolomics and measured AUC/CLr. Adding omics to the clinical model does not clearly improve prediction.','demo/#metformin','Voir le résultat','See the result',base=base)
        +'</div></section><section class="resource-band">'
        +'<div>'+t('Le projet est précisément intéressant parce que la réponse reste ouverte','The project matters precisely because the answer remains open','h2')+
        t('La plausibilité biologique existe, les contrats méthodologiques sont en place, mais la supériorité prédictive des omiques sur des références cliniques et PGx doit encore être démontrée sur des cohortes adaptées et indépendantes.',
        'Biological plausibility exists and methodological contracts are in place, but predictive superiority of omics over clinical and PGx references still has to be demonstrated in suitable independent cohorts.','p')+
        '</div>'+a('guide/','Voir comment le tester','See how to test it',base=base)+'</section>'
    )

def demo(base):
    metformin=REPO+'/blob/main/docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md'
    epa=REPO+'/blob/main/docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md'
    jarmusch=REPO+'/blob/main/docs/JARMUSCH14_PUBLIC_REAL_CLINICAL_PK_PGX_20261018.md'
    return ('<section class="small-hero"><p class="eyebrow">Démo scientifique / Scientific demo</p>'
        +t('Pourquoi Omics2PK ?', 'Why Omics2PK?','h1')
        +t('La question n’est pas « peut-on entraîner un modèle avec des omiques ? », mais « un état moléculaire mesuré avant l’exposition apporte-t-il une information PK utile au-delà de ce que l’on savait déjà ? »',
           'The question is not “can we train a model using omics?” but “does molecular state measured before exposure add useful PK information beyond what was already known?”','p','lede')+'</section>'
        '<section class="scientific-question" data-testid="central-hypothesis">'
        +t('Hypothèse centrale','Central hypothesis','p','eyebrow')
        +t('Chez les mêmes individus, une représentation omique pré-dose améliore-t-elle la prédiction tenue à l’écart d’une PK réellement mesurée par rapport à une référence clinique ± pharmacogénétique ?',
           'In the same people, does a predose omics representation improve held-out prediction of genuinely measured PK compared with a clinical ± pharmacogenetic reference?','h2')
        +t('Le mot important est « améliore ». Une association moléculaire ou une belle séparation de groupes ne suffit pas.',
           'The important word is “improve”. A molecular association or visually separated groups are not enough.','p')+'</section>'
        '<section class="evidence-step" id="epa"><div class="step-number">1</div><div>'
        +t('Première raison : la biologie rend l’hypothèse plausible','First reason: biology makes the hypothesis plausible','h2')
        +t('Dans une source U.S. EPA authentique, 100 foies humains possèdent ARN qPCR et activité enzymatique mesurée chez les mêmes donneurs. Après correction conservatrice pour huit CYP, trois comparaisons restent favorables à l’ARN : CYP1A2, CYP2B6 et activité CYP3A4.',
           'In an authentic U.S. EPA source, 100 human livers have matched qPCR RNA and measured enzyme activity. After conservative correction across eight CYP endpoints, three RNA comparisons remain favorable: CYP1A2, CYP2B6 and CYP3A4 activity.','p')
        +'<table><thead><tr><th>CYP</th><th>ΔMAE RNA − age/sex</th><th>Bonferroni 99.375%</th></tr></thead>'
        '<tbody><tr><td>CYP1A2</td><td>−5.589</td><td>[−10.310, −1.087]</td></tr>'
        '<tr><td>CYP2B6</td><td>−11.757</td><td>[−20.195, −3.224]</td></tr>'
        '<tr><td>CYP3A4 activity</td><td>−337.278</td><td>[−646.051, −37.048]</td></tr></tbody></table>'
        +'<div class="takeaway"><strong>'+t('Ce que cela montre','What this shows')+'</strong>'
        +t('Un état moléculaire peut contenir de l’information fonctionnelle. Cela ne prouve ni une AUC, ni une clairance systémique, ni une meilleure dose chez un patient.','Molecular state can contain functional information. This proves neither AUC, systemic clearance nor a better patient dose.','span')+'</div>'
        +a(epa,'Rapport source et multiplicité ↗','Source and multiplicity report ↗',base=base,external=True)
        +'</div></section>'
        '<section class="evidence-step" id="jarmusch"><div class="step-number">2</div><div>'
        +t('Deuxième raison : le potentiel génétique et l’état fonctionnel ne sont pas la même chose','Second reason: genetic potential and functional state are not the same thing','h2')
        +t('Dans 13 volontaires avec PK CYP3A répétée, une mesure fonctionnelle antérieure par sonde midazolam est numériquement meilleure que le génotype CYP disponible pour prédire la PK ultérieure : MAE 3,09465 contre 4,19763 sur l’échelle numérique source. L’intervalle descriptif du contraste est favorable, mais n13 reste très petit.',
           'In 13 volunteers with repeated CYP3A PK, a prior midazolam functional probe is numerically better than available CYP genotype for predicting later PK: MAE 3.09465 versus 4.19763 on the source numerical scale. The descriptive contrast interval is favorable, but n13 remains very small.','p')
        +'<div class="takeaway"><strong>'+t('Pourquoi c’est utile','Why this matters')+'</strong>'
        +t('Cela illustre le concept de phénotype dynamique que les omiques pourraient approcher. Mais une PK de sonde administrée est déjà un test fonctionnel pharmacologique : ce n’est pas une preuve que les omiques battent le génotype.','It illustrates the dynamic phenotype concept that omics may approximate. But an administered probe PK is already a pharmacological functional test; it is not evidence that omics outperform genotype.','span')+'</div>'
        +a(jarmusch,'Rapport clinique/PGx complet ↗','Full clinical/PGx report ↗',base=base,external=True)
        +'</div></section>'
        '<section class="evidence-step" id="metformin"><div class="step-number">3</div><div>'
        +t('Troisième raison : quand on teste enfin Omics→PK chez de vrais patients, le résultat peut être négatif','Third reason: when Omics→PK is finally tested in real patients, the result may be negative','h2')
        +t('Le benchmark metformine relie 32 personnes à une métabolomique pré-dose authentique de 384 variables et à deux critères PK mesurés. Les PCA, normalisations et modèles sont réappris dans chaque pli laissé de côté. Le test primaire compare clinique + omiques contre clinique seule.',
           'The metformin benchmark links 32 people to authentic 384-feature predose metabolomics and two measured PK endpoints. PCA, normalization and models are refit inside each held-out fold. The primary test compares clinical + omics with clinical alone.','p')
        +'<table><thead><tr><th>PK mesurée</th><th>Clinique MAE</th><th>Clinique + omiques MAE</th><th>ΔMAE</th></tr></thead><tbody>'
        '<tr><td>AUClast 0–24</td><td>2826.0472</td><td>2825.0240</td><td>−1.0232 (~0.04%)</td></tr>'
        '<tr><td>CLr</td><td>115.1113</td><td>118.1449</td><td>+3.0336 (worse)</td></tr></tbody></table>'
        +t('Les unités physiques absolues de la table source doivent encore être vérifiées ; l’interprétation porte donc sur la comparaison des erreurs, pas sur leur unité clinique. Pour l’AUC, l’intervalle descriptif de ΔMAE traverse largement zéro ; pour CLr, l’ajout omique est numériquement moins bon.',
           'Absolute physical units in the source table still require verification; interpretation therefore concerns relative prediction errors rather than clinical units. For AUC, the descriptive ΔMAE interval widely crosses zero; for CLr, adding omics is numerically worse.','p','note')
        +'<div class="takeaway"><strong>'+t('Résultat scientifique','Scientific result')+'</strong>'
        +t('La valeur ajoutée moléculaire n’est pas confirmée dans cette cohorte. C’est un résultat utile : Omics2PK doit comparer équitablement les modèles, pas produire artificiellement une histoire positive.','Molecular added value is not confirmed in this cohort. That is useful evidence: Omics2PK must compare models fairly rather than manufacture a positive story.','span')+'</div>'
        +a(metformin,'Benchmark humain original ↗','Original human benchmark ↗',base=base,external=True)
        +'</div></section>'
        '<section class="conclusion-panel" id="conclusion">'
        +t('Alors, pourquoi continuer ?', 'So why continue?','h2')
        +t('Parce que les trois étages ne donnent pas la même réponse : la plausibilité biologique existe ; un état fonctionnel peut apporter une information distincte du génotype ; mais le premier benchmark clinique omique du dépôt ne démontre pas encore de gain incrémental. La bonne question de recherche reste donc ouverte.',
           'Because the three evidence levels do not give the same answer: biological plausibility exists; functional state can add information distinct from genotype; but the repository’s first clinical omics benchmark does not yet demonstrate incremental gain. The research question therefore remains open.','p')
        +'<div class="decision-flow"><span>'+t('Plausible biologiquement','Biologically plausible')+'</span><b>→</b><span>'+t('Testable proprement','Can be tested rigorously')+'</span><b>→</b><span>'+t('Pas encore démontré en général','Not yet generally established')+'</span></div>'
        +'<div class="cta">'+a('guide/','Voir comment construire le test','See how to build the test',base=base)
        +'<a class="btn btn-outline" href="'+base+'analyze/" data-fr="État du futur outil Analyze" data-en="Status of the future Analyze tool">État du futur outil Analyze</a></div>'
        +'</section>')

def guide(base):
    return ('<section class="small-hero"><p class="eyebrow">Cours guidé / Guided course</p>'
        +t('Comment construire une étude Omics→PK correcte ?', 'How do you build a valid Omics→PK study?','h1')
        +t('Ce parcours suit la décision scientifique dans l’ordre. L’objectif n’est pas de choisir immédiatement un algorithme, mais d’empêcher les erreurs d’identité, de temporalité, de comparaison et de validation qui rendraient le résultat impossible à interpréter.',
           'This course follows the scientific decision in order. The goal is not to choose an algorithm immediately, but to prevent identity, timing, comparison and validation errors that would make the result uninterpretable.','p','lede')+'</section>'
        '<section class="course-map" aria-label="Course sequence">'
        '<span>1 · Question</span><b>→</b><span>2 · PK</span><b>→</b><span>3 · Temps</span><b>→</b><span>4 · Identité</span><b>→</b><span>5 · Référence</span><b>→</b><span>6 · Omique</span><b>→</b><span>7 · Validation</span><b>→</b><span>8 · Interprétation</span></section>'
        '<section class="course-step"><div class="step-number">1</div><div>'
        +t('Formuler la question avant de regarder les omiques','Define the question before looking at omics','h2')
        +t('Question recommandée : « chez des individus tenus à l’écart, l’état moléculaire mesuré avant l’exposition améliore-t-il la prédiction d’une PK réellement mesurée au-delà d’un modèle clinique ± PGx ? »',
           'Recommended question: “in held-out people, does molecular state measured before exposure improve prediction of genuinely measured PK beyond a clinical ± PGx model?”','p')
        +'<div class="do-dont"><div><strong>'+t('Correct','Correct')+'</strong>'+t('Comparaison incrémentale sur les mêmes personnes et le même critère PK.','Incremental comparison on the same people and same PK endpoint.','p')+'</div><div><strong>'+t('À éviter','Avoid')+'</strong>'+t('Chercher d’abord les variables qui corrèlent le mieux à la PK puis appeler cela validation.','Selecting the variables most correlated with PK first and then calling that validation.','p')+'</div></div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">2</div><div>'
        +t('Choisir une vraie cible PK','Choose a genuine PK target','h2')
        +t('Une preuve clinique Omics→PK exige une observation pharmacocinétique du médicament : AUC, CL/F, clairance, Cmax ou concentrations C(t) avec dose et temps. Une activité CYP, un métabolite formé in vitro ou un score de voie sont des preuves mécanistiques intermédiaires, pas de la PK clinique.',
           'Clinical Omics→PK evidence requires a pharmacokinetic observation of the drug: AUC, CL/F, clearance, Cmax, or C(t) concentrations with dose and time. CYP activity, in-vitro metabolite formation, or a pathway score are intermediate mechanistic evidence, not clinical PK.','p')
        +'<div class="method-note"><strong>'+t('Décision','Decision')+'</strong>'+t('Critère scalaire → comparaison prédictive simple possible. Concentrations–temps → passer à un modèle longitudinal/NLME distinct.','Scalar endpoint → a simpler predictive benchmark is possible. Concentration–time data → use a distinct longitudinal/NLME model.','span')+'</div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">3</div><div>'
        +t('Verrouiller la temporalité','Lock the temporal order','h2')
        +t('L’omique utilisée comme prédicteur doit être mesurée avant la dose ou avant le landmark PK que l’on veut prédire. Une mesure post-dose peut déjà contenir l’effet du médicament : l’utiliser comme « covariable pré-dose » crée une fuite d’information biologique.',
           'Omics used as a predictor must be measured before the dose or before the PK landmark being predicted. A post-dose assay may already contain drug effects: using it as a “predose covariate” creates biological information leakage.','p')
        +'<div class="timeline"><span>'+t('Omique pré-dose','Predose omics')+'</span><b>→</b><span>'+t('Dose index','Index dose')+'</span><b>→</b><span>'+t('PK observée','Observed PK')+'</span></div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">4</div><div>'
        +t('Prouver que ce sont les mêmes personnes','Prove that they are the same people','h2')
        +t('Le croisement doit utiliser un identifiant biologique authentifié ou un crosswalk source vérifié. Le tri des lignes, l’ordre des fichiers ou des labels ressemblants ne sont jamais une preuve d’identité. Les réplicats techniques restent attachés à la même personne.',
           'Joining datasets requires an authenticated biological identifier or source-verified crosswalk. Row order, file order or similar labels are never proof of identity. Technical replicates remain attached to the same person.','p')
        +'</div></section>'
        '<section class="course-step"><div class="step-number">5</div><div>'
        +t('Définir la référence que les omiques doivent battre','Define the baseline omics must beat','h2')
        +t('Au minimum, utilisez les covariables cliniques pertinentes réellement disponibles. Ajoutez la pharmacogénétique seulement si elle a été véritablement mesurée. Une absence de génotype doit être déclarée ; elle ne doit jamais être remplacée par un génotype simulé ou inféré sans validation.',
           'At minimum, use genuinely available relevant clinical covariates. Add pharmacogenetics only when it was actually measured. Missing genotype must be declared; it must never be replaced by simulated or unvalidated inferred genotype.','p')
        +'<div class="decision-grid"><article><strong>Baseline A</strong>'+t('Clinique','Clinical','p')+'</article><article><strong>Baseline B</strong>'+t('Clinique + PGx authentique','Clinical + authentic PGx','p')+'</article><article><strong>Test</strong>'+t('Référence + représentation omique','Baseline + omics representation','p')+'</article></div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">6</div><div>'
        +t('Choisir la représentation omique sans regarder le test','Choose the omics representation without peeking at the test set','h2')
        +t('Options raisonnables : variables directes régularisées, PCA apprise uniquement sur le train, score biologique fixé à l’avance, réseau/SWEET préspécifié, ou fusion multi-omique. Toute normalisation, sélection de variables, PCA, tuning ou apprentissage de score doit être réestimé dans le train de chaque pli.',
           'Reasonable options include regularized direct features, training-only PCA, a prespecified biological score, a prespecified network/SWEET representation, or multi-omics fusion. Normalization, feature selection, PCA, tuning or learned scores must be refit inside each training fold.','p')
        +'<div class="representation-row"><span>Direct</span><span>PCA</span><span>'+t('Score fixé','Fixed score')+'</span><span>SWEET / réseau</span><span>Multi-omics</span></div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">7</div><div>'
        +t('Valider au niveau de la personne, puis du site','Validate at the person level, then the site','h2')
        +t('Toutes les mesures d’un patient restent dans le même pli. Si plusieurs centres existent, une analyse leave-one-site-out teste mieux le transfert qu’un K-fold mélangeant les centres. Une cohorte externe réellement indépendante reste l’étape la plus forte.',
           'All measurements from one patient stay in the same fold. With multiple centers, leave-one-site-out better tests transportability than K-fold mixing centers. A genuinely independent external cohort remains the strongest step.','p')
        +'<div class="validation-ladder"><span>'+t('Patient tenu à l’écart','Held-out patient')+'</span><b>→</b><span>'+t('Centre tenu à l’écart','Held-out site')+'</span><b>→</b><span>'+t('Cohorte externe','External cohort')+'</span></div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">8</div><div>'
        +t('Interpréter le gain, y compris quand il n’existe pas','Interpret the gain, including when there is none','h2')
        +t('La sortie principale doit comparer les performances sur exactement les mêmes sujets : ΔMAE, ΔRMSE, biais, P20 et incertitude. Un modèle omique plus complexe qui n’améliore pas la référence est un résultat négatif informatif, pas un échec logiciel.',
           'The main output should compare performance on exactly the same people: ΔMAE, ΔRMSE, bias, P20 and uncertainty. A more complex omics model that does not improve the baseline is an informative negative result, not a software failure.','p')
        +'<div class="takeaway"><strong>'+t('Exemple réel du dépôt','Real repository example')+'</strong>'+t('Dans la cohorte metformine n32, clinique + métabolomique n’améliore pratiquement pas l’AUC et dégrade le MAE de CLr. La conclusion correcte est « valeur ajoutée non confirmée ».','In the n32 metformin cohort, clinical + metabolomics hardly changes AUC performance and worsens CLr MAE. The correct conclusion is “added value not confirmed”.','span')+'</div>'
        +'</div></section>'
        '<section class="course-step"><div class="step-number">9</div><div>'
        +t('Ne passer au PopPK qu’après avoir une question et des données compatibles','Move to PopPK only when the question and data support it','h2')
        +t('Si les données contiennent doses, temps et concentrations répétées, la covariable omique peut ensuite entrer dans un modèle NLME/PopPK avec IIV, diagnostics et comparaison de modèles. Une régression sur une AUC pré-estimée ne devient pas automatiquement un modèle PopPK.',
           'If the data contain doses, times and repeated concentrations, the omics covariate can then enter an NLME/PopPK model with IIV, diagnostics and model comparison. Regression on a pre-estimated AUC does not automatically become a PopPK model.','p')
        +'</div></section>'
        '<section class="course-step" id="checklist"><div class="step-number">10</div><div>'
        +t('Passer de MultiOmics à Omics2PK sans perdre la provenance','Handoff from MultiOmics to Omics2PK without losing provenance','h2')
        +t('MultiOmics prépare les matrices, identifiants et métadonnées moléculaires. Omics2PK ajoute ensuite, dans un environnement autorisé, les données cliniques/PGx/PK des mêmes personnes. Le contrat technique ne remplace jamais la vérification scientifique de la source.',
           'MultiOmics prepares molecular matrices, identifiers and metadata. Omics2PK then adds clinical/PGx/PK data for the same people in an authorized environment. The technical contract never replaces scientific source verification.','p')
        +'<ul class="checklist"><li>'+t('Même personne et même visite vérifiées','Same person and visit verified')+'</li><li>'+t('Omique antérieure à l’index PK','Omics precede the PK index')+'</li><li>'+t('PK réellement mesurée et unité documentée','PK genuinely measured and unit documented')+'</li><li>'+t('Référence clinique/PGx définie avant comparaison','Clinical/PGx baseline defined before comparison')+'</li><li>'+t('Prétraitement et sélection confinés au train','Preprocessing and selection confined to training')+'</li><li>'+t('Validation patient/site/externe explicitement nommée','Patient/site/external validation explicitly named')+'</li></ul>'
        +'<div class="cta">'+a('demo/','Revoir pourquoi ces garde-fous comptent','Review why these safeguards matter',base=base)
        +'<a class="btn btn-outline" href="'+MULTIOMICS+'" target="_blank" rel="noopener noreferrer" data-fr="Préparer les omiques dans MultiOmics ↗" data-en="Prepare omics in MultiOmics ↗">Préparer les omiques dans MultiOmics ↗</a></div>'
        +'</div></section>')

def learn(base):
    # Legacy route kept for old links; the navigation now points to /guide/.
    return guide(base)

def analyze(base):
    return ('<section class="small-hero"><p class="eyebrow">Analyze · en cours / in development</p>'
        +t('L’outil d’analyse Omics2PK n’est pas encore publié','The Omics2PK analysis tool is not published yet','h1')
        +t('Le dépôt contient déjà des moteurs, contrats, tests et benchmarks de recherche. Mais les exposer comme un outil web prêt à analyser des données réelles serait prématuré : le workflow public ne doit ni accepter des données patients, ni donner l’impression qu’une validation Omics→PK générale est acquise.',
           'The repository already contains research engines, contracts, tests and benchmarks. But exposing them as a web tool ready for real data would be premature: the public workflow must neither accept patient data nor imply that general Omics→PK validation has been achieved.','p','lede')+'</section>'
        '<section class="construction" data-testid="analyze-in-progress">'
        '<div class="status-pill">'+t('En cours de développement','In development')+'</div>'
        +t('Ce que devra faire Analyze','What Analyze will eventually do','h2')
        +'<div class="cards">'
        +card('01','Importer un échange MultiOmics qualifié','Import a qualified MultiOmics exchange','Conserver identités, couches, échelles, chronologie et provenance sans recopier les données cliniques protégées dans le site public.','Preserve identities, layers, scales, timing and provenance without copying protected clinical data into the public site.','guide/#checklist','Voir le contrat','See the contract',base=base)
        +card('02','Ajouter la référence PK dans un environnement autorisé','Add PK reference data in an authorized environment','Relier uniquement les mêmes personnes à des covariables cliniques, PGx et PK réellement mesurées.','Link only the same people to clinical covariates, PGx and genuinely measured PK in an authorized environment.','guide/','Voir la méthode','See the method',base=base)
        +card('03','Comparer équitablement les modèles','Compare models fairly','Rapporter le gain incrémental des représentations omiques contre les références sur les mêmes partitions patient/site.','Report incremental performance of omics representations against baselines on identical patient/site partitions.','demo/','Voir les preuves actuelles','See current evidence',base=base)
        +'</div>'
        +'<div class="callout">'+t('Pour l’instant, aucune zone de dépôt de données et aucun bouton « lancer l’analyse » ne sont affichés. Les utilisateurs de recherche avancés peuvent consulter le code et les protocoles reproductibles dans le dépôt GitHub.',
             'For now, no data-upload area and no “run analysis” button are shown. Advanced research users can inspect the reproducible code and protocols in the GitHub repository.','p')+'</div>'
        +'<div class="cta">'+a('demo/','Démo scientifique','Scientific demo',base=base)
        +'<a class="btn btn-outline" href="'+base+'guide/" data-fr="Cours guidé" data-en="Guided course">Cours guidé</a>'
        +'<a class="btn btn-outline" href="'+REPO+'" target="_blank" rel="noopener noreferrer" data-fr="Code de recherche ↗" data-en="Research code ↗">Code de recherche ↗</a></div>'
        +'</section>')

def disciplines(base):
    return ('<section class="small-hero"><p class="eyebrow">Ecosystem / Écosystème</p>'
         +t('Apprendre. Analyser. Intégrer.','Learn. Analyze. Integrate.','h1')
         +t('PMx Explain, MultiOmics et Omics2PK sont trois outils séparés reliés par une même logique de reproductibilité, de provenance et de niveaux de preuve explicites.',
           'PMx Explain, MultiOmics and Omics2PK are separate tools connected by shared principles of reproducibility, provenance and explicit evidence levels.','p','lede')+'</section>'
         +'<figure class="ecosystem-figure"><img src="'+base+'assets/ecosystem-graphical-abstract.svg" alt="PMx Explain for learning, MultiOmics for analysis and Omics2PK for PK integration">'
         +t('Parcours commun ; dépôts et responsabilités scientifiques séparés.','Shared journey; separate repositories and scientific responsibilities.','figcaption')+'</figure>'
         +'<div class="cards section">'
         +card('01 · Learn / Apprendre','PMx Explain','PMx Explain',
           'Comprendre PK/PD, PopPK, covariables et MIPD avant d’interpréter un signal biologique.',
           'Understand PK/PD, PopPK, covariates and MIPD before interpreting a biological signal.',
           PMX,'Ouvrir','Open',base=base)
         +card('02 · Analyze / Analyser','MultiOmics','MultiOmics',
           'Contrôler la qualité et produire des analyses, signatures et exports omiques traçables.',
           'Quality-control omics data and produce traceable analyses, signatures and exports.',
           MULTIOMICS,'Ouvrir','Open',base=base)
         +card('03 · Integrate / Intégrer','Omics2PK','Omics2PK',
           'Tester si des états moléculaires individuels améliorent une question PK au-delà des données cliniques et pharmacogénétiques.',
           'Test whether individual molecular states improve a PK question beyond clinical and pharmacogenetic information.',
           'demo/','Voir la recherche','See the research',base=base)
         +'</div>'
         +'<section class="ecosystem-boundary">'
         +t('La chaîne scientifique','Scientific handoff','p','eyebrow')
         +t('Échange moléculaire vérifié ≠ bénéfice PK démontré','Verified molecular exchange ≠ demonstrated PK benefit','h2')
         +t('MultiOmics peut fournir des mesures et métadonnées individuelles traçables. Omics2PK exige ensuite les mêmes sujets avec chronologie, covariables de référence, doses et concentrations authentiques. La compatibilité technique seule ne valide ni une covariable ni une recommandation de dose.',
           'MultiOmics can provide traceable individual measurements and metadata. Omics2PK then requires the same subjects with timing, reference covariates and genuine dose/concentration observations. Technical compatibility alone validates neither a covariate nor a dosing recommendation.','p')
         +'</section>')

def resources(base):
    refs=[
        ('Clinical Omics→PK: metformin n32 predose metabolomics + measured AUC/CLr',REPO+'/blob/main/docs/E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md'),
        ('Clinical functional phenotype vs PGx: Jarmusch CYP probes',REPO+'/blob/main/docs/JARMUSCH14_PUBLIC_REAL_CLINICAL_PK_PGX_20261018.md'),
        ('100 human livers – original U.S. EPA paired RNA/activity','https://doi.org/10.23719/1519197'),
        ('Original EPA source and multiplicity-corrected numerical results',REPO+'/blob/main/docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md'),
        ('Fair single- vs multi-omics → observed PK protocol',REPO+'/blob/main/docs/MULTIMODAL_PK_ABLATION_SOFTWARE_QUALIFICATION_20261019.md'),
        ('Native MultiOmics → Omics2PK exchange contract',REPO+'/blob/main/docs/NATIVE_MULTIOMICS_INTEGRATION.md'),
        ('Research reproducibility and PK protocol',REPO+'/blob/main/docs/E2_NCA_E3_POPPK_PBPK_E4_VALIDATION_PROTOCOL_20261016.md'),
        ('Mavoglurant longitudinal NLME pilot and its limits',REPO+'/blob/main/docs/MAVOGLURANT_2C_NLME_IIV_CL_PILOT.md'),
        ('Source code / license',REPO),
    ]
    return ('<section class="small-hero"><p class="eyebrow">Resources</p>'
      +t('Sources, protocoles et reproductibilité','Sources, protocols and reproducibility','h1')
      +t('Chaque démonstration renvoie à sa source et à son statut de preuve. Aucun fichier individuel protégé n’est redistribué.',
         'Every demonstration links back to its source and evidence status. No protected individual records are redistributed.','p','lede')+'</section><section class="section"><ol class="references">'
      +''.join('<li><a href="'+escape(url,quote=True)+'" target="_blank" rel="noopener noreferrer">'+escape(name)+' ↗</a></li>' for name,url in refs)
      +'</ol></section>')

def document(page,body,base):
    titles={'home':'Omics2PK – omiques vers pharmacocinétique',
            'demo':'Démo scientifique – Omics2PK',
            'guide':'Cours guidé – Omics2PK',
            'learn':'Cours guidé – Omics2PK',
            'analyze':'Analyze – en cours – Omics2PK',
            'ecosystem':'Ecosystem – Omics2PK','disciplines':'Ecosystem – Omics2PK','resources':'Resources – Omics2PK'}
    return ('<!doctype html><html lang="fr"><head>'
       '<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
       '<meta name="theme-color" content="#fbfaf8"><meta name="referrer" content="strict-origin-when-cross-origin">'
       '<meta name="description" content="Omics2PK : démonstration scientifique et méthodologie guidée pour tester rigoureusement les liens entre omiques et pharmacocinétique.">'
       '<title>'+titles[page]+'</title><link rel="stylesheet" href="'+base+'assets/styles.css">'
       '<script defer src="'+base+'assets/app.js"></script></head><body>'
       '<a href="#content" class="skip-link" data-fr="Aller au contenu" data-en="Skip to content">Aller au contenu</a>'
       '<div class="app">'+nav(base,page)+'<main id="content" class="page-shell">'+body+'</main>'
       '<footer class="site-footer"><div class="footer-shell"><div><strong>Omics2PK</strong>'
       +t('Tester rigoureusement si l’état moléculaire améliore la compréhension ou la prédiction PK.','Rigorously test whether molecular state improves PK understanding or prediction.','p')
       +'</div><nav aria-label="Liens de pied de page">'+a(PMX,'PMx Explain','PMx Explain',base=base,external=True)
       +a(MULTIOMICS,'MultiOmics','MultiOmics',base=base,external=True)
       +a(REPO,'Code source ↗','Source code ↗',base=base,external=True)
       +'</nav></div></footer><div class="disclaimer" role="note">'
       +t('Plateforme de recherche exploratoire — aucun outil clinique de dose.',
          'Exploratory research platform — no clinical dosing tool.','span')
       +'</div></div></body></html>')

def render(out):
    out=Path(out)
    mapping={'home':(Path('index.html'),'.',home),
             'demo':(Path('demo/index.html'),'..',demo),
             'guide':(Path('guide/index.html'),'..',guide),
             'learn':(Path('learn/index.html'),'..',learn),
             'analyze':(Path('analyze/index.html'),'..',analyze),
             'ecosystem':(Path('ecosystem/index.html'),'..',disciplines),
             'disciplines':(Path('disciplines/index.html'),'..',disciplines),
             'resources':(Path('resources/index.html'),'..',resources)}
    for page,(path,depth,func) in mapping.items():
        base='./' if depth=='.' else '../'
        dest=out/path
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(document(page,func(base),base),encoding='utf-8')
    asset=Path(__file__).resolve().parent/'assets'
    import shutil
    shutil.copytree(asset,out/'assets',dirs_exist_ok=True)
    (out/'.nojekyll').write_text('',encoding='utf-8')
    (out/'robots.txt').write_text('User-agent: *\nAllow: /\n',encoding='utf-8')
    canonical=('','demo','guide','analyze','ecosystem','disciplines','resources')
    (out/'sitemap.xml').write_text(
      '<?xml version="1.0" encoding="UTF-8"?>'
      '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
      +''.join('<url><loc>'+ROOT+(page+'/' if page else '')+'</loc></url>' for page in canonical)
      +'</urlset>',encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='site/dist')
    render(parser.parse_args().out)

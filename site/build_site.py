"""Render the Omics2PK static research portal, same editorial design as PMx Explain.

No JavaScript framework, tracking, external fonts, data upload, or clinical claims.
Static project Pages paths: /Omics2PK/{learn,analyze,disciplines,resources}/
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

def nav(base,page):
    groups = [
      ("learn","Learn","Learn",[
         ("learn/","Apprendre les méthodes","Learn the methods",False),
         ("learn/#sources","Comprendre les preuves","Understand evidence",False)
      ]),
      ("analyze","Analyze","Analyze",[
         ("analyze/","Explorer les résultats validés","Explore validated results",False),
         ("analyze/#methodes","Méthodes reproductibles","Reproducible methods",False)
      ]),
      ("disciplines","Other Disciplines","Other Disciplines",[
         (PMX,"PMx Explain ↗","PMx Explain ↗",True),
         (MULTIOMICS,"MultiOmics ↗","MultiOmics ↗",True),
         ("disciplines/","Écosystème des disciplines","Discipline ecosystem",False)
      ]),
      ("resources","Ressources","Resources",[
         ("resources/","Références et documentation","References and documentation",False),
         (REPO,"Code source et rapports ↗","Source code and reports ↗",True),
         (MAIN,"Site principal ↗","Main website ↗",True)
      ])
    ]
    details=[]
    icons={"learn":"▤","analyze":"◫","disciplines":"◇","resources":"▦"}
    for key,fr,en,links in groups:
        inner="".join(a(path,lfr,len_,base=base,external=ext) for path,lfr,len_,ext in links)
        details.append('<details class="goal" name="site-objective" data-testid="goal-'+key+'">'
          '<summary'+(' class="active"' if page==key else '')+'><span class="nav-icon" aria-hidden="true">'+icons[key]+'</span>'
          +t(fr,en)+'<span class="chevron" aria-hidden="true">⌄</span></summary><div class="goal-links">'+inner+'</div></details>')
    return ('<header class="site-header"><a class="logo" href="'+base+'" aria-label="Omics2PK, accueil" data-testid="logo-link">'
        '<span class="mark" aria-hidden="true">oP</span><span class="word">Omics<em>2PK</em></span></a>'
        '<button class="burger" type="button" aria-label="Ouvrir le menu" aria-controls="site-nav" aria-expanded="false" data-testid="nav-toggle"><span></span><span></span><span></span></button>'
        '<nav id="site-nav" aria-label="Navigation principale" data-testid="primary-nav">'
        +''.join(details)
        +'<div class="tools"><button type="button" id="lang-switch" aria-label="Switch language" aria-pressed="false">FR / EN</button>'
        '<button type="button" id="theme-switch" aria-label="Changer le thème" aria-pressed="false">◐ <span data-fr="Thème" data-en="Theme">Thème</span></button></div></nav></header>')

def card(kicker,title_fr,title_en,desc_fr,desc_en,href,link_fr,link_en,*,base):
    return ('<article class="card"><div class="card-kicker">'+escape(kicker)+'</div>'
        +t(title_fr,title_en,'h3')+t(desc_fr,desc_en,'p')
        +'<div class="card-link">'+a(href,link_fr+' →',link_en+' →',base=base,external=href.startswith('http'))+'</div></article>')

def home(base):
    return (
        '<section class="hero"><p class="eyebrow">Omics2PK · research laboratory</p>'
        +t("Relier les omiques à la pharmacocinétique.",
           "Connecting omics to pharmacokinetics.",'h1')
        +t("Des signatures moléculaires mesurées aux modèles de métabolisme et, à terme, à la pharmacocinétique de population. Une plateforme de recherche transparente : jeux de données originaux, analyses reproductibles et limites explicites.",
           "From measured molecular signatures to drug metabolism and, ultimately, population pharmacokinetics. An open research platform with original datasets, reproducible analyses and explicit limitations.",'p','lede')
        +'<div class="cta">'+a("analyze/","Explorer les validations","Explore the evidence",base=base)
        +'<a class="btn btn-outline" href="'+base+'learn/" data-fr="Comprendre la méthode" data-en="Learn the method">Comprendre la méthode</a></div>'
        +t("Recherche exploratoire — aucun modèle de dose clinique validé.", "Exploratory research only — no clinically validated dosing model.",'p','hero-note')+'</section>'
        +'<section class="objectives" aria-labelledby="objective-title">'+t("Quel est votre objectif ?", "What is your objective?",'h2')+
        '<div class="objectives-grid">'+
        a('learn/','Comprendre les liens biologiques','Understand biological relationships',base=base)+
        a('analyze/','Examiner des résultats réels','Inspect measured results',base=base)+
        a('resources/','Reproduire une analyse','Reproduce an analysis',base=base)+
        a('disciplines/','Explorer les outils associés','Explore related tools',base=base)+
        '</div></section><section class="section"><p class="eyebrow">Evidence first</p>'
        +t('Preuves disponibles et limites', 'What is established — and what is not','h2')
        +'<div class="cards">'
        +card('E2 · verified','100 foies humains et huit CYP','100 human livers, eight CYPs','Expression qPCR et activité enzymatique mesurées sur les mêmes donneurs. Trois signaux restent favorables après correction pour huit comparaisons.','Matched donor qPCR and measured enzyme activity. Three favorable signals after correction across eight outcomes.','analyze/#epa','Lire les résultats','Read the results',base=base)
        +card('E2 · verified','Deux inhibiteurs de tyrosine kinase','Two tyrosine-kinase inhibitors','Tucatinib et imatinib : mêmes 21 donneurs microsomaux. Prédiction de la clairance in vitro à partir des protéines CYP ; ce ne sont pas deux cohortes indépendantes.','Tucatinib and imatinib: the same 21 microsomal donors. In-vitro clearance prediction from CYP proteins; not two independent cohorts.','analyze/#tki','Voir la comparaison','See the comparison',base=base)
        +card('E2 · measured','Sulfaméthoxazole, NAT1 et NAT2','Sulfamethoxazole, NAT1 and NAT2','26 hépatocytes-donneurs originaux, modèles sur n26/n21/n20 et comparaison exacte sur 15 donneurs. Aucun gain statistiquement démontré.','26 original hepatocyte donors, n26/n21/n20 models and identical-subject 15-donor comparison. No statistically established improvement.','analyze/#smx','Lire les résultats','Read measured results',base=base)
        +'</div></section><section class="resource-band">'
        +'<div>'+t('De la preuve enzymatique à la PK clinique','From enzyme evidence to clinical PK','h2')+
        t('La prochaine étape exige des données omiques, des doses, des horaires et des concentrations réellement associées chez les mêmes personnes.',
        'The next step needs genuinely linked same-patient omics, dosing history, timestamps and measured drug concentrations.','p')+
        '</div>'+a('resources/','Consulter les jeux de données','Review source datasets',base=base)+'</section>'
    )

def learn(base):
    return ('<section class="small-hero"><p class="eyebrow">Learn / Apprendre</p>'
        +t('Des omiques à la PK : ce que mesure chaque étape','From omics to PK: what each stage measures','h1')
        +t('La validation biologique ne remplace pas la validation chez des patients exposés au médicament.',
           'Biological validation does not substitute for validation in treated patients.','p','lede')+'</section>'
        '<section class="section" id="sources"><div class="cards">'
        +card('01','Préparer une matrice omique','Prepare an omics matrix','Identifier les donneurs, distinguer qPCR, protéomique et métabolomique, vérifier les unités, l’origine et la chronologie.','Link original donors; identify qPCR, proteomics and metabolomics, units, provenance and timing.',MULTIOMICS,'Découvrir MultiOmics','Explore MultiOmics',base=base)
        +card('02','Relier omiques et métabolisme','Connect omics and metabolism','Comparer ARN/protéines et activité CYP/NAT ou CLint au même donneur, avec validation donneur exclu.','Match CYP/NAT activity or CLint to molecular features in the same biological donor, using donor-held-out validation.','analyze/','Voir les preuves','View evidence',base=base)
        +card('03','Passer aux observations PK','Advance to observed PK','Exiger dose, temps, formulation, concentrations et covariables individuelles ; comparer aux modèles NLME/popPK pertinents.','Require observed doses, time, formulation, concentrations and individual covariates, then compare to appropriate NLME/PopPK baselines.','resources/','Lire le protocole','Read the protocol',base=base)
        +'</div></section><section class="callout">'+t('Important : la formation d’un métabolite in vitro n’est ni une concentration plasmatique, ni une AUC, ni une recommandation de dose.',
          'Important: in-vitro metabolite formation is neither a plasma concentration, an AUC, nor a dose recommendation.','p')+'</section>')

def analyze(base):
    source1='https://github.com/rberrah/Omics2PK/blob/main/docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md'
    source2='https://github.com/rberrah/Omics2PK/blob/main/docs/MEASURED_TKI_TWO_DRUG_ORIGINAL_RESULTS_20261017.md'
    issue='https://github.com/rberrah/Omics2PK/issues/38'
    return ('<section class="small-hero"><p class="eyebrow">Analyze / Analyser</p>'
        +t('Résultats réels, scripts et limites','Measured results, source scripts and limitations','h1')
        +t('Cette page présente des analyses pré-calculées et vérifiées ; elle n’est pas un moteur en ligne de calcul PK ou de prédiction thérapeutique.',
           'This page shows precomputed, source-verified research results; it is not an online PK analysis engine or clinical dosing predictor.','p','lede')+'</section>'
        '<section class="section evidence" id="epa"><p class="eyebrow">Source original verified · n = 100</p>'
        +t('ARN hépatique → activité CYP','Liver RNA → CYP activity','h2')
        +t('Après correction conservatrice de Bonferroni sur huit CYP, les intervalles exploratoires excluent zéro pour trois enzymes.',
           'After conservative Bonferroni correction across eight CYP endpoints, exploratory intervals exclude zero for three enzymes.','p')
        +'<table><thead><tr><th>CYP</th><th>ΔMAE RNA − age/sex</th><th>Bonferroni 99.375%</th></tr></thead>'
        '<tbody><tr><td>CYP1A2</td><td>−5.589</td><td>[−10.310, −1.087]</td></tr>'
        '<tr><td>CYP2B6</td><td>−11.757</td><td>[−20.195, −3.224]</td></tr>'
        '<tr><td>CYP3A4 activity</td><td>−337.278</td><td>[−646.051, −37.048]</td></tr></tbody></table>'
        +t('Unité : pmol/min/mg de protéines. Quatre comparaisons nominalement favorables au seuil 95 %, trois après correction ; CYP2E1 est ponctuellement moins bon avec ARN.',
           'Units: pmol/min/mg protein. Four nominal 95% intervals but only three after correction; CYP2E1 is directionally worse with RNA.','p','note')
        +a(source1,'Méthodes, chiffres et CI source ↗','Methods, results and source CI ↗',base=base,external=True)
        +'</section><section class="section evidence" id="tki"><p class="eyebrow">Source original verified · same 21 donors</p>'
        +t('Protéines CYP → clairance microsomale','CYP proteins → microsomal clearance','h2')
        +'<table><thead><tr><th>Outcome</th><th>n</th><th>Baseline MAE</th><th>Protein MAE</th></tr></thead><tbody>'
        '<tr><td>Tucatinib CLint</td><td>21</td><td>6.246</td><td>2.041</td></tr>'
        '<tr><td>Imatinib quantified CLint</td><td>17</td><td>28.870</td><td>14.166</td></tr></tbody></table>'
        +t('Clairances microsomales exprimées en µL/min/mg. Les deux études utilisent les mêmes donneurs ; quatre clairances imatinib ne sont pas quantifiables.',
           'Microsomal clearance in µL/min/mg. Both studies share the same donors; four imatinib clearances are nonquantifiable.','p','note')
        +a(source2,'Rapport de qualification ↗','Qualification report ↗',base=base,external=True)
        +'</section><section class="section evidence" id="smx"><p class="eyebrow">SOURCE-GATED · no real numerical model result</p>'
        +t('NAT1/NAT2 → formation de NA-SMX','NAT1/NAT2 → NA-SMX formation','h2')
        +t('Recalcul depuis le Tableau 1 original : activités NAT (n21) MAE 0,03903 contre référence 0,04230 ; protéines NAT2 (n21) 0,04039. Comparaison stricte sur 15 mêmes donneurs : activités 0,03798, protéines 0,04382, génotype 0,03927 ; référence 0,04191. Tous les intervalles exploratoires des gains recouvrent zéro.',
           'Reanalysis of the original Table 1: NAT probe activities (n21) MAE 0.03903 vs baseline 0.04230; NAT2 proteins (n21) 0.04039. Strict identical-15-donor comparison: activities 0.03798, proteins 0.04382, genotype 0.03927; baseline 0.04191. All exploratory improvement intervals cross zero.','p')
        +t('Ces résultats proviennent de la transcription locale de la table de l’article, non d’un téléchargement XML automatisé en CI. Pas de données individuelles redistribuées.',
           'These findings derive from local transcription of the article source table, not from a CI-downloaded XML. No individual donor records are redistributed.','p','note')
        +a('https://github.com/rberrah/Omics2PK/blob/main/docs/SMX_NAT_26_ORIGINAL_MEASURED_LODO_RESULTS_20261017.md','Rapport des valeurs originales ↗','Original measured results report ↗',base=base,external=True)
        +' · '+a(issue,'Provenance et reproduction ↗','Source provenance and reproduction ↗',base=base,external=True)
        +'</section><section class="callout" id="methodes">'+t('Tous les chiffres affichés concernent la biologie et la PK in vitro, pas la prédiction clinique personnalisée.',
             'All displayed numbers concern in-vitro metabolism, not individualized clinical PK prediction.','p')+'</section>')

def disciplines(base):
    return ('<section class="small-hero"><p class="eyebrow">Other Disciplines</p>'
         +t('Un même écosystème, trois spécialisations','One ecosystem, three specialties','h1')
         +t('Une navigation harmonisée avec PMx Explain et MultiOmics ; chaque outil conserve sa mission et ses limites méthodologiques.',
           'Shared navigation and design with PMx Explain and MultiOmics; each tool retains its research scope and limitations.','p','lede')+'</section>'
         +'<div class="cards section">'
         +card('PK/PD','PMx Explain','PMx Explain','Cours, animations et outils de pharmacométrie ; point d’entrée pédagogique principal.','Pharmacometrics lessons, interactive visualizations and research tools.',PMX,'Visiter','Open site',base=base)
         +card('Omics','MultiOmics','MultiOmics','Analyse transcriptomique, protéomique et métabolomique : extraction reproductible des signatures moléculaires.','Transcriptomic, proteomic and metabolomic analyses to derive reproducible molecular signatures.',MULTIOMICS,'Visiter','Open site',base=base)
         +card('Molecular PK','Omics2PK','Omics2PK','Validation du lien omiques → métabolisme → PK ; analyses de recherche, non utilisation clinique.','Evidence connecting omics to metabolism and PK; research-only, not a clinical tool.','analyze/','Voir la recherche','See the research',base=base)
         +'</div>')

def resources(base):
    refs=[
        ('100 human livers – original U.S. EPA paired RNA/activity','https://doi.org/10.23719/1519197'),
        ('Original EPA source and multiplicity-corrected numerical results',REPO+'/blob/main/docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md'),
        ('Tucatinib HLM donor original article','https://doi.org/10.1016/j.dmd.2025.100061'),
        ('Imatinib HLM donor original article','https://doi.org/10.1016/j.dmd.2025.100196'),
        ('Sulfamethoxazole NAT1/NAT2 PHH source','https://doi.org/10.1016/j.dmd.2025.100195'),
        ('Original 26-donor SMX numerical source results',REPO+'/blob/main/docs/SMX_NAT_26_ORIGINAL_MEASURED_LODO_RESULTS_20261017.md'),
        ('Sulfamethoxazole original-data issue #38',REPO+'/issues/38'),
        ('Research reproducibility and PK protocol',REPO+'/blob/main/docs/E2_NCA_E3_POPPK_PBPK_E4_VALIDATION_PROTOCOL_20261016.md'),
        ('Source code / license',REPO),
    ]
    return ('<section class="small-hero"><p class="eyebrow">Resources</p>'
      +t('Sources, protocoles et reproductibilité','Sources, protocols and reproducibility','h1')
      +t('Chaque résultat renvoie à la publication, au dépôt original ou à la validation source. Aucun fichier individuel protégé n’est redistribué.',
         'Every result links to the article, original repository or source qualification. No protected individual donor records are redistributed.','p','lede')+'</section><section class="section"><ol class="references">'
      +''.join('<li><a href="'+escape(url,quote=True)+'" target="_blank" rel="noopener noreferrer">'+escape(name)+' ↗</a></li>' for name,url in refs)
      +'</ol></section>')

def document(page,body,base):
    titles={'home':'Omics2PK – omiques vers pharmacocinétique',
            'learn':'Learn – Omics2PK','analyze':'Analyze – Omics2PK',
            'disciplines':'Other Disciplines – Omics2PK','resources':'Resources – Omics2PK'}
    return ('<!doctype html><html lang="fr"><head>'
       '<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
       '<meta name="theme-color" content="#fbfaf8"><meta name="referrer" content="strict-origin-when-cross-origin">'
       '<meta name="description" content="Omics2PK : plateforme de recherche reproductible reliant omiques, métabolisme et pharmacocinétique.">'
       '<title>'+titles[page]+'</title><link rel="stylesheet" href="'+base+'assets/styles.css">'
       '<script defer src="'+base+'assets/app.js"></script></head><body>'
       '<a href="#content" class="skip-link" data-fr="Aller au contenu" data-en="Skip to content">Aller au contenu</a>'
       '<div class="app">'+nav(base,page)+'<main id="content" class="page-shell">'+body+'</main>'
       '<footer class="site-footer"><div class="footer-shell"><div><strong>Omics2PK</strong>'
       +t('Reproductibilité scientifique des ponts omiques–PK.','Reproducible science linking omics and PK.','p')
       +'</div><nav aria-label="Liens de pied de page">'+a(PMX,'PMx Explain','PMx Explain',base=base,external=True)
       +a(MULTIOMICS,'MultiOmics','MultiOmics',base=base,external=True)
       +a(REPO,'Code source ↗','Source code ↗',base=base,external=True)
       +'</nav></div></footer><div class="disclaimer" role="note">'
       +t('Outil de recherche exploratoire — aucune validation clinique de dose.',
          'Exploratory research tool — no clinically validated dosing.','span')
       +'</div></div></body></html>')

def render(out):
    out=Path(out)
    mapping={'home':(Path('index.html'),'.',home),
             'learn':(Path('learn/index.html'),'..',learn),
             'analyze':(Path('analyze/index.html'),'..',analyze),
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
    (out/'sitemap.xml').write_text(
      '<?xml version="1.0" encoding="UTF-8"?>'
      '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
      +''.join('<url><loc>'+ROOT+(page+'/' if page else '')+'</loc></url>' for page in ('','learn','analyze','disciplines','resources'))
      +'</urlset>',encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='site/dist')
    render(parser.parse_args().out)

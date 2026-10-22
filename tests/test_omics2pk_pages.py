"""Static site contract: scientific demo, guided course, in-progress Analyze page."""
import importlib.util
import re
import tempfile
import unittest
from pathlib import Path


class Omics2PKPagesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/"site"/"build_site.py"
        spec=importlib.util.spec_from_file_location("omics2pk_site_build",path)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        cls.generated=tempfile.TemporaryDirectory()
        cls.root=Path(cls.generated.name)
        module.render(cls.root)

    @classmethod
    def tearDownClass(cls):
        cls.generated.cleanup()

    def pages(self):
        return ["index.html","demo/index.html","guide/index.html","learn/index.html",
                "analyze/index.html","ecosystem/index.html","disciplines/index.html",
                "resources/index.html"]

    def test_static_pages_declare_their_public_canonical_url(self):
        origin="https://rberrah.github.io/Omics2PK/"
        for path in self.pages():
            html=(self.root/path).read_text(encoding="utf-8")
            route=path.replace("index.html","")
            if route=="learn/":
                route="guide/"
            expected=origin+route
            self.assertEqual(html.count('rel="canonical"'),1,path)
            self.assertIn(f'<link rel="canonical" href="{expected}">',html,path)

    def test_pages_keep_shared_navigation_and_learn_subpages(self):
        for page in self.pages():
            html=(self.root/page).read_text(encoding="utf-8")
            self.assertEqual(html.count('<details class="goal"'),4,page)
            for goal in ("learn","analyze","disciplines","resources"):
                self.assertIn('data-testid="goal-'+goal+'"',html,page)
            self.assertIn('data-fr="Apprendre"',html,page)
            self.assertIn('data-fr="Analyze"',html,page)
            self.assertIn('data-fr="Écosystème"',html,page)
            self.assertIn('data-fr="Ressources"',html,page)
            self.assertIn('data-fr="Objectif"',html,page)
            self.assertIn('data-fr="Démo guidée"',html,page)
            self.assertIn("site-header",html)
            self.assertIn('data-testid="logo-link"',html)
            self.assertIn('data-testid="language-toggle"',html)
            self.assertIn('data-testid="theme-toggle"',html)
            self.assertIn('data-testid="nav-toggle"',html)
            self.assertIn('rel="noopener noreferrer"',html)
            self.assertIn('id="content"',html)
            self.assertIn('data-fr=',html)
            self.assertIn('data-en=',html)
            self.assertNotIn("<form",html.lower())

    def test_local_routes_really_resolve_from_each_page(self):
        for page in self.pages():
            html=(self.root/page).read_text(encoding="utf-8")
            for href in re.findall(r'href="([^"]+)"',html):
                if href.startswith(('http:', 'https:', '#')):
                    continue
                clean=href.split("#")[0].split("?")[0]
                if not clean:
                    continue
                target=(self.root/Path(page).parent/clean).resolve()
                self.assertTrue(target.is_relative_to(self.root.resolve()), (page,href))
                self.assertTrue(target.is_file() or (target.is_dir() and (target/"index.html").exists()),
                                (page,href))
            for path in re.findall(r'src="([^"]+)"',html):
                target=(self.root/Path(page).parent/path).resolve()
                self.assertTrue(target.exists(), (page,path))

    def test_scientific_demo_answers_why_without_overclaiming(self):
        demo=(self.root/"demo/index.html").read_text(encoding="utf-8")
        self.assertIn('data-testid="central-hypothesis"',demo)
        for term in ("−5.589","−11.757","−337.278",
                     "−0.2724","−0.4721","−0.0905","−0.0084",
                     "3,09465","4,19763",
                     "7.343413","9.888019","4.601193","4.748841",
                     "2826.0472","2825.0240","115.1113","118.1449"):
            self.assertIn(term,demo)
        self.assertIn("valeur ajoutée moléculaire n’est pas confirmée",demo)
        self.assertIn("Omics2PK conserve ce résultat négatif",demo)
        self.assertIn('id="wegler"',demo)
        self.assertIn("Benchmark Wegler donneur-protéine→CLint",demo)
        self.assertIn("elle n’est pas un substitut universel",demo)
        self.assertIn("Pourquoi Omics2PK",demo)
        self.assertIn("Apprendre · Objectif",demo)
        self.assertNotIn("lancer l’analyse",demo.lower())

    def test_guided_course_teaches_how_in_scientific_order(self):
        guide=(self.root/"guide/index.html").read_text(encoding="utf-8")
        legacy=(self.root/"learn/index.html").read_text(encoding="utf-8")
        for term in ("Comment construire une étude Omics→PK correcte",
                     "Verrouiller la temporalité",
                     "Prouver que ce sont les mêmes personnes",
                     "Définir la référence que les omiques doivent battre",
                     "Valider au niveau de la personne, puis du site",
                     "Ne passer au PopPK"):
            self.assertIn(term,guide)
        self.assertIn('id="checklist"',guide)
        self.assertIn("Même personne et même visite vérifiées",guide)
        self.assertIn("Apprendre · Démo guidée / Learn · Guided demo",guide)
        self.assertIn("Apprendre · Démo guidée / Learn · Guided demo",legacy)

    def test_analyze_is_explicitly_in_progress_and_has_no_analyzer(self):
        analyze=(self.root/"analyze/index.html").read_text(encoding="utf-8")
        self.assertIn('data-testid="analyze-in-progress"',analyze)
        self.assertIn("En cours de développement",analyze)
        self.assertIn("n’est pas encore publié",analyze)
        self.assertIn("aucune zone de dépôt de données",analyze)
        self.assertNotIn("<form",analyze.lower())
        self.assertNotIn('type="file"',analyze.lower())
        self.assertNotIn('type="submit"',analyze.lower())

    def test_ecosystem_destinations_and_graphical_abstract(self):
        for page in ("index.html","ecosystem/index.html","disciplines/index.html"):
            html=(self.root/page).read_text(encoding="utf-8")
            self.assertIn('https://rberrah.github.io/',html)
            self.assertIn('https://rberrah.github.io/multiomics/',html)
            self.assertIn('https://github.com/rberrah/Omics2PK',html)
        ecosystem=(self.root/"ecosystem/index.html").read_text(encoding="utf-8")
        self.assertIn('Ecosystem / Écosystème',ecosystem)
        self.assertIn('PMx Explain',ecosystem)
        self.assertIn('MultiOmics',ecosystem)
        self.assertIn('Omics2PK',ecosystem)
        self.assertIn('ecosystem-graphical-abstract.svg',ecosystem)
        self.assertIn('../demo/',ecosystem)
        self.assertTrue((self.root/"assets/ecosystem-graphical-abstract.svg").is_file())

    def test_theme_css_and_behavior_and_no_tracking(self):
        css=(self.root/"assets/styles.css").read_text()
        js=(self.root/"assets/app.js").read_text()
        for token in ("#fbfaf8","#14161a","#a4441f","#e8845c",
                      "59px","@media(max-width:1100px)",
                      ".evidence-step",".course-step",".status-pill"):
            self.assertIn(token,css)
        self.assertIn("prefers-reduced-motion",css)
        self.assertIn("localStorage",js)
        self.assertIn("Escape",js)
        for name in ("inter-latin-wght-normal.woff2", "source-serif-4-latin-wght-normal.woff2", "jetbrains-mono-latin-wght-normal.woff2"):
            self.assertIn(name,css)
        self.assertIn("pharmacometrie-language",js)
        self.assertIn("pk-theme",js)
        self.assertIn("modes = ['system', 'light', 'dark']",js)
        self.assertNotIn("fetch(",js)
        self.assertNotIn("google-analytics",css+js)

    def test_built_site_sitemap_uses_canonical_new_routes(self):
        self.assertTrue((self.root/".nojekyll").exists())
        sitemap=(self.root/"sitemap.xml").read_text()
        self.assertEqual(sitemap.count("<loc>"),7)
        self.assertIn("https://rberrah.github.io/Omics2PK/demo/",sitemap)
        self.assertIn("https://rberrah.github.io/Omics2PK/guide/",sitemap)
        self.assertIn("https://rberrah.github.io/Omics2PK/analyze/",sitemap)
        self.assertIn("https://rberrah.github.io/Omics2PK/ecosystem/",sitemap)
        self.assertNotIn("https://rberrah.github.io/Omics2PK/learn/",sitemap)


if __name__ == "__main__":
    unittest.main()

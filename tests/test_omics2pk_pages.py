"""Static site contract: multi-page links, research evidence, nav/theme language."""
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

    def test_five_pages_with_identical_four_goal_dropdowns(self):
        pages=["index.html","learn/index.html","analyze/index.html",
               "disciplines/index.html","resources/index.html"]
        for page in pages:
            html=(self.root/page).read_text(encoding="utf-8")
            self.assertEqual(html.count('<details class="goal"'),4,page)
            for goal in ("learn","analyze","disciplines","resources"):
                self.assertIn('data-testid="goal-'+goal+'"',html,page)
            self.assertIn("site-header",html)
            self.assertIn('data-testid="logo-link"',html)
            self.assertIn('id="lang-switch"',html)
            self.assertIn('id="theme-switch"',html)
            self.assertIn('data-testid="nav-toggle"',html)
            self.assertIn('rel="noopener noreferrer"',html)
            self.assertIn('id="content"',html)
            self.assertIn('data-fr=',html)
            self.assertIn('data-en=',html)
            self.assertNotIn("<form",html.lower())

    def test_local_routes_really_resolve_from_each_page(self):
        for page in ["index.html","learn/index.html","analyze/index.html",
                     "disciplines/index.html","resources/index.html"]:
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

    def test_original_source_claims_and_explicit_abstention(self):
        evidence=(self.root/"analyze/index.html").read_text(encoding="utf-8")
        for term in ("−5.589","−11.757","−337.278","6.246","2.041",
                     "28.870","14.166","SOURCE-GATED"):
            self.assertIn(term,evidence)
        self.assertIn("original donor records have not yet generated a validated MAE",evidence)
        self.assertIn("no clinically validated dosing", evidence.lower())

    def test_external_discipline_destinations(self):
        for page in ("index.html","disciplines/index.html"):
            html=(self.root/page).read_text(encoding="utf-8")
            self.assertIn('https://rberrah.github.io/',html)
            self.assertIn('https://rberrah.github.io/multiomics/',html)
            self.assertIn('https://github.com/rberrah/Omics2PK',html)

    def test_theme_css_and_behavior_and_no_tracking(self):
        css=(self.root/"assets/styles.css").read_text()
        js=(self.root/"assets/app.js").read_text()
        for token in ("#fbfaf8","#14161a","#a4441f","#e8845c",
                      "59px","@media(max-width:1100px)"):
            self.assertIn(token,css)
        self.assertIn("prefers-reduced-motion",css)
        self.assertIn("localStorage",js)
        self.assertIn("Escape",js)
        self.assertNotIn("fetch(",js)
        self.assertNotIn("google-analytics",css+js)

    def test_built_site_sitemap_and_404_free_project_routes(self):
        self.assertTrue((self.root/".nojekyll").exists())
        sitemap=(self.root/"sitemap.xml").read_text()
        self.assertEqual(sitemap.count("<loc>"),5)
        self.assertIn("https://rberrah.github.io/Omics2PK/analyze/",sitemap)


if __name__ == "__main__":
    unittest.main()

"""Ensure the public research graphical abstract is versioned, accessible and safe."""
from pathlib import Path
import re
import unittest
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SVG = ROOT / "docs/assets/omics2pk-graphical-abstract.svg"
READ = ROOT / "README.md"


class GraphicalAbstractReleaseTests(unittest.TestCase):
    def test_svg_is_valid_accessible_editable_local_vector(self):
        root = ET.parse(SVG).getroot()
        ns = "{http://www.w3.org/2000/svg}"
        self.assertEqual(root.tag, ns + "svg")
        self.assertEqual(root.attrib.get("role"), "img")
        self.assertEqual(root.attrib.get("viewBox"), "0 0 1500 800")
        self.assertEqual(root.attrib.get("aria-labelledby"), "ga-title ga-desc")
        self.assertTrue((root.find(ns + "title").text or "").startswith("Omics2PK"))
        self.assertIn("not yet demonstrated", root.find(ns + "desc").text)
        self.assertGreater(len(root.findall(".//" + ns + "text")), 12)
        for el in root.iter():
            self.assertNotIn(el.tag.rsplit("}", 1)[-1].lower(),
                             {"script", "foreignobject", "image"})
            for attr, value in el.attrib.items():
                self.assertFalse(attr.endswith("href") and
                                 str(value).startswith(("http://", "https://")))
        self.assertLess(SVG.stat().st_size, 100_000)

    def test_readme_has_correct_explanatory_status_and_local_graphic(self):
        readme = READ.read_text(encoding="utf-8")
        self.assertIn("](docs/assets/omics2pk-graphical-abstract.svg)", readme)
        self.assertIn("no independently validated gain in tacrolimus", readme.lower())
        self.assertIn("clinical + authentic CYP3A4/5 genotype", readme)
        self.assertIn("not clinical software", readme.lower())
        self.assertTrue((ROOT / "docs/history/README_CHRONOLOGICAL_ARCHIVE_20261019.md").is_file())
        for target in re.findall(r"\]\((docs/[^)]+)\)", readme):
            self.assertTrue((ROOT / target.split("#")[0]).is_file(), target)


if __name__ == "__main__":
    unittest.main()

"""Crawler of the public chemical database PubChem."""
from typing import List

from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException, NoSuchElementException
import pandas as pd

from scripts.crawlers.abstract_crawler import AbstractCrawler


class PubchemCrawler(AbstractCrawler):
    """
    Class for the Pubchem crawler. The class uses Selenium for information extraction.
    The crawler uses webdriver for Chrome. The class Follows the :class:`AbstractCrawler` interface.
    """
    @classmethod
    def get_classname(cls) -> str:
        """
        Get the name of the class.
        :return: name of the class as a string.
        """
        return cls.__name__

    COLUMNS = ["Phototoxic", "Name", "Smiles", "Molecular Weight", "Hydrogen Bond Donor Count",
               "Hydrogen Bond Acceptor Count", "Rotatable Bond Count", "Exact Mass",
               "Monoisotopic Mass", "Topological Polar Surface Area", "Heavy Atom Count",
               "Formal Charge", "Complexity", "Isotope Atom Count",
               "Defined Atom Stereocenter Count", "Undefined Atom Stereocenter Count",
               "Defined Bond Stereocenter Count", "Undefined Bond Stereocenter Count",
               "Covalently-Bonded Unit Count", "Compound Is Canonicalized"]

    @classmethod
    def __extract_properties(cls, driver: webdriver.Chrome) -> List[list]:
        """
        On the webpage for a chemical on Pubchem, find a table containing computed properties
        and extract information from this table.
        :param driver: Chrome webdriver.
        :return: list of lists containing Chem+Phys properties.
        """
        properties = WebDriverWait(driver, 60).until(EC.presence_of_element_located(
            (By.ID, "Computed-Properties")))
        table = properties.find_element_by_tag_name("table")
        property_values = []
        for row in table.find_elements_by_xpath(".//tr"):
            new = [td.text for td in row.find_elements_by_xpath(".//td")]
            if len(new) != 0:
                property_values.append(new)
        return property_values

    @classmethod
    def crawl(cls, prefs: dict) -> None:
        """
        Extracts relevant information from the Pubchem database. Results are saved into a .csv file.
        :param prefs: preferences for the Chrome webdriver.
        :return: None
        """
        options = webdriver.ChromeOptions()
        options.add_experimental_option("prefs", prefs)
        driver = webdriver.Chrome("crawlers/dependencies/chromedriver.exe", options=options)

        smiles_df = pd.DataFrame(columns=cls.COLUMNS)

        df = pd.read_csv("crawlers/cache/chem.csv")

        # Querying pubchem for every chemical from the paper
        for name, phototoxic in zip(df['name'], df['phototoxic']):
            driver.get("https://pubchem.ncbi.nlm.nih.gov/")

            # Input and submit chem name
            driver.find_element_by_xpath(
                "//*[@class='flex-grow-1 p-xsm-left p-xsm-top p-xsm-bottom b-right b-light']").\
                find_element_by_tag_name("input").send_keys(name)
            driver.find_element_by_xpath("//*[@class='button width-2em height-2em lh-1']").click()

            # Wait for chem page to load
            compound = WebDriverWait(driver, 60).until(EC.presence_of_element_located(
                (By.XPATH, "//*[@class='f-medium p-sm-top p-sm-bottom f-1125']")))
            compound.click()

            # Extract Chem+Phys properties from the table on the page
            property_values = cls.__extract_properties(driver)
            cas_code = driver.find_element_by_id("CAS").find_element_by_tag_name("p").text

            # Attempt to find direct link to CAS to extract SMILES code
            try:
                driver.find_element_by_xpath("//*[@data-label='Reference: CAS Open']").click()
                driver.find_element_by_xpath(
                    "//a[contains(@href,'/source/CAS Common Chemistry')]//following::a[1]").click()
                smiles_code = WebDriverWait(driver, 5).until(EC.presence_of_element_located(
                    (By.XPATH, "//*[@class='smile']//following::p[1]"))).text
            # If link isn't provided, query CAS manually
            except (TimeoutException, NoSuchElementException):
                driver.get("https://commonchemistry.cas.org/")
                WebDriverWait(driver, 60).until(EC.presence_of_element_located(
                    (By.XPATH, ("//*[@class='search-input query-input ng-untouched "
                                "ng-pristine ng-valid']")))).send_keys(cas_code)
                driver.find_element_by_xpath("//*[@class='btn btn-primary btn-search']").click()
                try:
                    WebDriverWait(driver, 5).until(EC.presence_of_element_located(
                        (By.XPATH, "//*[@class='result-content'][1]"))).click()
                    smiles_code = WebDriverWait(driver, 5).until(EC.presence_of_element_located(
                        (By.XPATH, "//*[@class='smile']//following::p[1]"))).text
                except TimeoutException:
                    continue

            new_row = {'Phototoxic': phototoxic, 'Name': name, 'Smiles': smiles_code}
            for value in property_values:
                # Remove Angstrom unit symbol
                if value[0] == "Topological Polar Surface Area":
                    value[1] = value[1].split(' ')[0]
                new_row[value[0]] = value[1]

            smiles_df = smiles_df.append(new_row, ignore_index=True)

        smiles_df.to_csv("data/crawler_output/pubchem.csv")
        driver.close()

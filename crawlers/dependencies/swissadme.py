from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import os
import time
import pandas as pd


def is_file_downloaded(filename, timeout=5):
    end_time = time.time() + timeout
    while not os.path.exists(filename):
        time.sleep(1)
        if time.time() > end_time:
            print("File not found within time")
            return False

    if os.path.exists(filename):
        print("File found")
        return True


def crawl(prefs):
    download_path = prefs['download.default_directory'] + r'\swissadme.csv'

    options = webdriver.ChromeOptions()

    options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome('crawlers/dependencies/chromedriver.exe', options=options)
    driver.set_page_load_timeout(1200)

    pubchem_df = pd.read_csv('data/crawler_output/pubchem.csv')
    smiles_list = pubchem_df['Smiles']
    name_list = pubchem_df['Name']

    driver.get("http://www.swissadme.ch/index.php")
    text_area = WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.ID, "smiles")))

    for smiles in smiles_list:
        text_area.send_keys(smiles + '\n')

    try:
        driver.find_element_by_id('submitButton').click()
    except Exception:
        pass

    WebDriverWait(driver, 1200).until(EC.presence_of_element_located((By.XPATH, f"//*[@name='{len(smiles_list)}']")))
    driver.find_element_by_xpath("//a[contains(@href, 'swissadme.csv')]").click()

    if is_file_downloaded(download_path):
        adme_df = pd.read_csv(download_path)
        adme_df = adme_df.rename(columns={'Molecule': 'Name', 'Canonical SMILES': 'Smiles'})
        adme_df['Name'] = name_list
        adme_df['Smiles'] = smiles_list
        os.remove(download_path)
    else:
        driver.close()
        print(f'Swissadme download failed')
        return

    adme_df.to_csv('data/crawler_output/swiss.csv')
    driver.close()

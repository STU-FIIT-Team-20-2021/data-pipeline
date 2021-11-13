from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException, NoSuchElementException
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
    # Specify your own download folder
    download_path = prefs['download.default_directory'] + r'\SwissTargetPrediction.csv'

    options = webdriver.ChromeOptions()

    options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome('crawlers/dependencies/chromedriver.exe', options=options)

    predictions_df = pd.DataFrame(columns=['Compound', 'Target', 'Common name', 'Uniprot ID', 'ChEMBL ID', 'Target Class', 'Probability*', 'Known actives (3D/2D)'])

    df = pd.read_csv('data/crawler_output/pubchem.csv')
    smiles_list = df['Smiles']
    name_list = df['Name']

    for name, smiles in zip(name_list, smiles_list):
        driver.get('http://www.swisstargetprediction.ch/')

        try:
            smiles = smiles.split('.')[0]
        except:
            pass

        try:
            driver.find_element_by_id('smilesBox').send_keys(smiles)
            driver.find_element_by_id('submitButton').click()
        except NoSuchElementException:
            continue

        try:
            elem = WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.XPATH, "//*[@class='dt-button buttons-csv buttons-html5']")))
            elem.click()
        except TimeoutException:
            continue
        except NoSuchElementException:
            continue


        if is_file_downloaded(download_path):
            predictions_df = pd.concat([predictions_df, pd.read_csv(download_path)], axis=0, ignore_index=True)
            predictions_df['Compound'] = predictions_df['Compound'].fillna(name)
            os.remove(download_path)
        else:
            print(f'{smiles} targets download failed')

    predictions_df.to_csv('data/crawler_output/target.csv')
    driver.close()








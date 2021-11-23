"""Defines an abstract class for crawlers."""
import os
import time
from abc import ABC, abstractmethod


class AbstractCrawler(ABC):
    """Abstract class that is used as an interface for crawler classes."""
    @classmethod
    def is_file_downloaded(cls, filename: str, timeout: int = 5) -> bool:
        """
        Checks whether the required file is present.
        :param filename: path to the file.
        :param timeout: time until file search stops.
        :return: boolean stating whether the file was found or not.
        """
        end_time = time.time() + timeout
        while not os.path.exists(filename):
            time.sleep(1)
            if time.time() > end_time:
                print("File not found within time")
                return False

        if os.path.exists(filename):
            print("File found")
            return True

    @classmethod
    @abstractmethod
    def crawl(cls, prefs: dict) -> None:
        """
        Abstract method, has to be implemented. It should contain the core functionality
        of the crawler.
        :param prefs: preferences for the webdriver.
        :return: None
        """

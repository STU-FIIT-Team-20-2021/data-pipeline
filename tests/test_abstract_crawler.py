"""Unit tests for the AbstractCrawler class."""
import os
import threading

from scripts.crawlers.abstract_crawler import AbstractCrawler

TEST_FILENAME = "test_file.txt"


def create_file(f_name: str) -> None:
    """
    Utility function for file creation.
    :param f_name: file path
    :return: none
    """
    f = open(f_name, "x")
    f.close()


def test_file_downloaded():
    """
    Unit test. Checks AbstractCrawler's is_file_downloaded method.
    Scenario: file is present.
    Assert: method's result is true.
    """
    create_file(TEST_FILENAME)
    true_result = AbstractCrawler.is_file_downloaded(TEST_FILENAME)
    os.remove(TEST_FILENAME)
    assert true_result is True


def test_file_not_downloaded():
    """
    Unit test. Checks AbstractCrawler's is_file_downloaded method.
    Scenario: file is not present.
    Assert: method's result is false.
    """
    false_result = AbstractCrawler.is_file_downloaded(TEST_FILENAME)
    assert false_result is False


def test_file_downloaded_during_function():
    """
    Unit test. Checks AbstractCrawler's is_file_downloaded method.
    Scenario: file is created during the method's execution.
    Assert: method's result is true.
    """
    timer = threading.Timer(2, lambda: create_file(TEST_FILENAME))
    timer.start()
    true_result = AbstractCrawler.is_file_downloaded(TEST_FILENAME)
    os.remove(TEST_FILENAME)
    assert true_result is True

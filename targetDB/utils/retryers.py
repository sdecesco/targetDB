#!/usr/bin/env python

import urllib.request as urllib

import requests
import time
import pandas as pd


class NetworkError(RuntimeError):
	pass


def retryer(max_retries=10, timeout=5):
	def wraps(func):
		request_exceptions = (
			requests.exceptions.Timeout,
			requests.exceptions.ConnectionError,
			requests.exceptions.HTTPError, urllib.URLError, urllib.HTTPError
		)

		def inner(*args, **kwargs):
			for i in range(max_retries):
				try:
					result = func(*args, **kwargs)
				except request_exceptions:
					time.sleep(timeout * i)
					continue
				else:
					return result
			else:
				raise NetworkError

		return inner

	return wraps


def retryer_pubmed(max_retries=10, timeout=5):
	def wraps(func):
		request_exceptions = (
			requests.exceptions.Timeout,
			requests.exceptions.ConnectionError,
			requests.exceptions.HTTPError, urllib.URLError, urllib.HTTPError
		)

		def inner(*args, **kwargs):
			for i in range(max_retries):
				try:
					result = func(*args, **kwargs)
				except request_exceptions:
					time.sleep(timeout * i)
					continue
				else:
					return result
			else:
				return pd.DataFrame(data=None)

		return inner

	return wraps

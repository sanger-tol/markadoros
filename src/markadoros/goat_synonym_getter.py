import json

import requests
from loguru import logger


class GOATSynonymGetter:
    """
    Query GoaT (Genomes on a Tree) for a given taxon to retrieve a list of valid
    synonyms.
    """

    def __init__(self) -> None:
        self.goat_url = "https://goat.genomehubs.org/api/v2/lookup"

    def query_goat(self, search_term: str) -> requests.Response:
        headers = {"Accept": "application/json"}
        search_string = f"?searchTerm={search_term}&result=taxon&size=1&suggestSize=0"
        goat_response = requests.get(self.goat_url + search_string, headers=headers)

        return goat_response

    def process_goat_result(
        self, goat_response: requests.Response, taxon: str
    ) -> list[str]:
        try:
            goat_data = goat_response.json()

            if not goat_data.get("status").get("success"):
                logger.warning("GOAT lookup failed - no synonyms will be used!")
                return []

            first_result = goat_data.get("results")[0].get("result")

            if first_result.get("scientific_name") != taxon:
                logger.warning(
                    f"GOAT lookup returned a different taxon: {first_result.get('scientific_name')}."
                )
                synonyms = [
                    x.get("name")
                    for x in first_result.get("taxon_names")
                    if x.get("class") in ["scientific name", "synonym"]
                ]

                if taxon in synonyms:
                    synonyms.remove(taxon)
                    return synonyms
                else:
                    return []

            else:
                synonyms = [
                    x.get("name")
                    for x in first_result.get("taxon_names")
                    if x.get("class") == "synonym"
                ]

            logger.info(f"Found the following synonyms: {', '.join(synonyms)}")

            return synonyms

        except json.JSONDecodeError:
            logger.warning("Failed to decode GOAT response - no synonyms will be used!")
            return []

    def get_synonyms(self, taxon_name: str) -> list[str]:
        goat_response = self.query_goat(taxon_name)

        if goat_response.status_code != 200:
            logger.warning(
                f"GOAT lookup failed (status code '{str(goat_response.status_code)}') - no synonyms will be used!"
            )
            return []

        return self.process_goat_result(goat_response, taxon_name)

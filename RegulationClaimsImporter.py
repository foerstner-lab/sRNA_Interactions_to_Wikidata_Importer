import pywikibot


class RegulationClaimsImporter:
    def __init__(self, claims_df):
        self.claims_df = claims_df

    def import_claims(self):
        imports_count = 0
        for claim_index, claim in self.claims_df.iterrows():
            regulator_qid = claim['regulator_QID']
            property_pid = claim['property']
            target_qid = claim['target_QID']
            site = pywikibot.Site("wikidata", "wikidata")
            repo = site.data_repository()
            item = pywikibot.ItemPage(repo, regulator_qid)
            item.get()
            if item.claims:
                if property_pid in item.claims:
                    if item.claims[property_pid][0].getTarget() == target_qid:
                        continue
            claim = pywikibot.Claim(repo, property_pid)
            target = pywikibot.ItemPage(repo, target_qid)
            claim.setTarget(target)
            item.addClaim(claim, summary=u'Adding claim')
            imports_count += 1
        return imports_count

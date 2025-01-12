# OPTIMADE, MPContribs databases

!pip install optimade[http_client]
!pip install mpcontribs-client

from optimade.client import OptimadeClient

client = OptimadeClient(
    include_providers={"mp"},
)
# client.get('elements HAS "Bi"')
client.list_properties("structures")

from mpcontribs.client import Client

client = Client()
client.available_query_params()  # print list of available query parameters

query = {"formula__contains": "Bi", "project": "qmof"}
fields = ["id", "identifier", "formula", "structures"]
data = client.query_contributions(
    query=query, fields=fields, sort="", paginate=False
)

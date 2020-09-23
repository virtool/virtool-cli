import os
import time
import json
import asyncio
from Bio import Entrez
from concurrent.futures.thread import ThreadPoolExecutor

Entrez.email  = os.environ["NCBI_EMAIL"]
API_KEY = os.environ["NCBI_API_KEY"]


async def run(src_path: str):
    executor = ThreadPoolExecutor(max_workers=10)
    loop = asyncio.get_event_loop()

    paths = get_paths(src_path)

    coros = list()

    for path in paths[:5]:
        name = get_name_from_path(path)
        
        if name:
            coro = loop.run_in_executor(executor, fetch_taxid, name)
            coros.append(coro)

            await asyncio.sleep(0.3)
        else:
            paths.remove(path)

    results = await asyncio.gather(*coros)
    
    
    names = []
    taxids = []
    # parse through results
    for name, taxid in results:
        names.append(name)
        taxids.append(taxid)
    
    for name, taxid, path, in zip(names, taxids, paths):
        with open(os.path.join(path, "otu.json"), 'r+') as f:
            otu = json.load(f)

            if otu["name"] == name:
                otu["taxid"] = taxid
                f.seek(0)
                json.dump(otu, f, indent=4)
    
    if not results:
        print("All OTUs are up to date")

    print(f"Retrieved {len(taxids)} taxids for {len(names)} OTUs")
        


def fetch_taxid(name):
    handle = Entrez.esearch(db="taxonomy", term=name, api_key=API_KEY)
    record = Entrez.read(handle)

    try:
        taxid = record["IdList"][0]
        print(f"Received: {name} {taxid}")
    except IndexError:
        print(f"Could not receive: {name}")
        taxid = None


    return name, taxid


def get_paths(src_path: str):
    alpha_paths = os.listdir(src_path)
    paths = []

    for alpha in alpha_paths:
        if alpha == "meta.json":
            continue

        otu_paths = [os.path.join(src_path, alpha, otu) for otu in os.listdir(os.path.join(src_path, alpha))]
        
        for otu in otu_paths:
            paths.append(otu)

    return paths


def get_name_from_path(path):
    with open(os.path.join(path, "otu.json"), 'r') as f:
        otu = json.load(f)

        try:
            if otu["taxid"]: return None
        except KeyError:
            return otu["name"]


def taxid(src_path):
    loop = asyncio.get_event_loop()
    loop.run_until_complete(run(src_path))

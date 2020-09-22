import os
import config
import time
import json
import asyncio
import concurrent.futures
from Bio import Entrez
from concurrent.futures.thread import ThreadPoolExecutor


def taxid(src_path: str, force_update: bool):
    otu_paths = get_paths(src_path)
    Entrez.email = config.email

    start = time.time()
    with ThreadPoolExecutor(max_workers=10) as executor:
        executor.map(fetch, otu_paths, [force_update]*len(otu_paths))
        executor.shutdown(wait=True)
    
    print(str(time.time()-start))
    #asyncio.run(run(otu_paths, force_update))
    

# async def run(otu_paths: list, force_update: bool):
    
#     tasks = []
#     for path in otu_paths:
#         task = asyncio.create_task(fetch(path, force_update))
#         tasks.append(task)

#     await asyncio.gather(*tasks)
#     #await asyncio.wait(tasks)
#     #asyncio.run(tasks)


def fetch(otu_path: str, force_update: bool):

    with open(os.path.join(otu_path, "otu.json"), "r") as f:
        otu = json.load(f)
    
    print("Finding " + otu["name"])
    #await asyncio.sleep(0.2)
    handle = Entrez.esearch(db="taxonomy", term=otu["name"], api_key=config.API_KEY)
    record = Entrez.read(handle)
    print(otu['name'] + " done")

    if "taxid" not in otu or force_update:
        try:
            otu["taxid"] = record["IdList"][0]
        except IndexError:
            otu["taxid"] = None


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


if __name__ == "__main__":
    taxid("tests/files/src", True)

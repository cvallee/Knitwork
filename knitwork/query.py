import mrich
import logging
from mrich import print

import json
import time
import asyncio
from rdkit.Chem import MolFromSmiles
from neo4j import GraphDatabase, AsyncGraphDatabase, Query

from .config import CONFIG
from .tools import load_sig_factory, calc_pharm_fp


def check_config():
    graph_vars = ["GRAPH_LOCATION", "GRAPH_USERNAME", "GRAPH_PASSWORD"]
    missing = []
    for var in graph_vars:
        if var not in CONFIG:
            missing.append(var)
    if missing:
        mrich.error("Configuration missing:", missing)
        raise ValueError(f"Configuration missing: {missing}")


def get_driver():
    check_config()
    return GraphDatabase.driver(
        CONFIG["GRAPH_LOCATION"],
        auth=(CONFIG["GRAPH_USERNAME"], CONFIG["GRAPH_PASSWORD"]),
    )


async def aget_driver():
    check_config()
    return AsyncGraphDatabase.driver(
        CONFIG["GRAPH_LOCATION"],
        auth=(CONFIG["GRAPH_USERNAME"], CONFIG["GRAPH_PASSWORD"]),
    )


async def arun_query(query, timeout=None):
    driver = await aget_driver()
    async with driver:
        async with driver.session() as session:
            if timeout:
                result = await session.run(Query(query, timeout=timeout))
            else:
                result = await session.run(Query(query, timeout=timeout))
            records = [record async for record in result]
            return records


def run_query(query, timeout=None):
    driver = get_driver()
    with driver:
        with driver.session() as session:
            if timeout:
                result = session.run(Query(query), timeout=timeout)
            else:
                result = session.run(Query(query))
            records = [record for record in result]
            return records


async def aget_subnodes(
    smiles: str,
    num_hops: int = 25,
    terminal_nodes: bool = CONFIG["FRAGMENT_TERMINAL_SUBNODES"],
    prevent_cylces: bool = True,
    progress=None,
    task=None,
    timeout=None,
):
    """
    Get subnodes for a given node (retrieve using SMILES)

    :param smiles: SMILES string for node to retrieve subnodes
    :param terminal_subnodes: whether to only return 'terminal' subnodes (can't be broken down further)
    :return: list of unique subnode SMILES
    """

    if prevent_cylces:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})
        CALL apoc.path.expandConfig(a, {
          relationshipFilter: 'FRAG>',
          minLevel: 0,
          maxLevel: %(num_hops)d,
          uniqueness: 'RELATIONSHIP_GLOBAL'
        }) YIELD path
        WITH last(nodes(path)) AS f, relationships(path) AS e
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }
    else:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})-[e:FRAG*0..%(num_hops)d]->(f:F2)
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }

    if terminal_nodes:
        query = query + " WHERE NOT ()-[:FRAG]-(f)-[:FRAG]->()"

    query = query + " RETURN f"

    records = await arun_query(query, timeout=timeout)
    subnodes = [record["f"]["smiles"] for record in records]

    if progress:
        progress.update(task, advance=1)

    return set(subnodes)


async def aget_synthons(
    smiles: str,
    num_hops: int = 25,
    terminal_nodes: bool = CONFIG["FRAGMENT_TERMINAL_SYNTHONS"],
    prevent_cylces: bool = True,
    progress=None,
    task=None,
    timeout=None,
):
    """
    Get constituent synthons (compounds added or removed during transformation) for a given node SMILES.
    [Xe] denotes the attachment point.

    :param smiles: SMILES string of node to retrieve synthons
    :param terminal_synthons: whether to return 'terminal' synthons, i.e. can't be broken down more
    :return: list of constituent synthon SMILES strings
    """

    if prevent_cylces:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})
        CALL apoc.path.expandConfig(a, {
          relationshipFilter: 'FRAG>',
          minLevel: 0,
          maxLevel: %(num_hops)d,
          uniqueness: 'RELATIONSHIP_GLOBAL'
        }) YIELD path
        WITH last(nodes(path)) AS b, relationships(path) AS e
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }
    else:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})-[e:FRAG*0..%(num_hops)d]->(b:F2)
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }

    if terminal_nodes:
        query = query + " WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()"

    query = query + " RETURN e[-1] AS edge"

    records = await arun_query(query, timeout=timeout)
    edges = [edge for record in records if (edge := record["edge"])]

    synthons = set()

    for edge in edges:
        for p in ["prop_synthon", "prop_core"]:
            i = edge[p]
            if not i:
                continue

            if i.count("Xe") != 1:
                continue

            synthons.add(i)

    if progress:
        progress.update(task, advance=1)

    return synthons


async def aget_r_groups(
    smiles: str,
    num_hops: int = 25,
    prevent_cylces: bool = True,
    progress=None,
    task=None,
    timeout=None,
):

    if prevent_cylces:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})
        CALL apoc.path.expandConfig(a, {
          relationshipFilter: 'FRAG>',
          minLevel: 0,
          maxLevel: %(num_hops)d,
          uniqueness: 'RELATIONSHIP_GLOBAL'
        }) YIELD path
        WITH last(nodes(path)) AS b, relationships(path) AS e
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        AND e[-1].prop_synthon contains '[Xe]'
        AND NOT e[-1].prop_synthon=e[-2].prop_synthon
        RETURN e[-1].prop_synthon AS synthon, e[-2].prop_synthon AS r_group;
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }
    else:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})-[e:FRAG*0..%(num_hops)d]->(b:F2)
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        AND e[-1].prop_synthon contains '[Xe]'
        AND NOT e[-1].prop_synthon=e[-2].prop_synthon
        RETURN e[-1].prop_synthon AS synthon, e[-2].prop_synthon AS r_group;
        """ % {
            "smiles": smiles,
            "num_hops": num_hops,
        }

    records = await arun_query(query, timeout=timeout)
    results = []
    for record in records:
        results.append((record["synthon"], record["r_group"]))

    if progress:
        progress.update(task, advance=1)

    return results


def get_pure_expansions(
    smiles: str,
    synthon: str,
    num_hops: int = 2,
    limit: int = 5,
    index: int | None = None,
    prevent_cylces: bool = True,
    cache_dir=None,
    cached_only=False,
    timeout=None,
):

    if cache_dir:
        cache_file = cache_dir / f"pure_{smiles}_{synthon}_{num_hops}_{limit}.json"
        if cache_file.exists():
            logging.info(f"Using cache {index} {smiles} {synthon}")
            return json.load(open(cache_file, "rt"))
        elif cached_only:
            return None

    if prevent_cylces:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})
        CALL apoc.path.expandConfig(a, {
          relationshipFilter: '<FRAG',
          minLevel: 0,
          maxLevel: %(num_hops)d,
          uniqueness: 'RELATIONSHIP_GLOBAL'
        }) YIELD path
        WITH last(nodes(path)) AS b
        MATCH (b)<-[e:FRAG]-(c:Mol)
        WHERE e.prop_synthon = '%(synthon)s'
        RETURN c.smiles AS smi, c.cmpd_ids AS ids
        """ % {
            "smiles": smiles,
            "synthon": synthon,
            "num_hops": num_hops
            }
    else:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
        WHERE e.prop_synthon='%(synthon)s'
        WITH c.smiles AS smi, c.cmpd_ids AS ids
        RETURN smi, ids
        """ % {
            "smiles": smiles,
            "synthon": synthon,
            "num_hops": num_hops
        }

    if limit:
        query = query + f" LIMIT {limit}"

    logging.info(f"Starting pure expansion {index} {smiles} {synthon}")

    try:
        records = run_query(query, timeout=timeout)
    except Exception as e:
        mrich.error(index, e)
        raise Exception(f"{smiles=} {synthon=} {e}")

    results = []
    for record in records:
        results.append((record["ids"], record["smi"]))

    if cache_dir:
        json.dump(results, open(cache_file, "wt"), indent=2)

    logging.info(f"Success {index} {smiles} {synthon} #results: {len(results)}")

    return results


def get_impure_expansions(
    smiles: str,
    synthon: str,
    num_hops: int = 2,
    limit: int = 5,
    index: int | None = None,
    prevent_cylces: bool = True,
    cache_dir=None,
    cached_only=False,
    timeout=None,
):

    if cache_dir:
        cache_file = cache_dir / f"impure_{smiles}_{synthon}_{num_hops}_{limit}.json"
        if cache_file.exists():
            logging.info(f"Using cache {index} {smiles} {synthon}")
            return json.load(open(cache_file, "rt"))
        elif cached_only:
            return None

    logging.info(f"Starting impure expansion {index} {smiles} {synthon}")

    sig_factory = load_sig_factory(
        fdef_file=CONFIG["FINGERPRINT_FDEF"],
        max_point_count=CONFIG["FINGERPRINT_MAXPOINTCOUNT"],
        bins=json.loads(CONFIG["FINGERPRINT_BINS"]),
    )

    vector = calc_pharm_fp(MolFromSmiles(synthon), sig_factory, as_str=True)

    threshold = CONFIG["KNITWORK_SIMILARITY_THRESHOLD"]
    metric = CONFIG["KNITWORK_SIMILARITY_METRIC"]

    if prevent_cylces:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})
        CALL apoc.path.expandConfig(a, {
          relationshipFilter: '<FRAG',
          minLevel: 0,
          maxLevel: %(num_hops)d,
          uniqueness: 'RELATIONSHIP_GLOBAL'
        }) YIELD path
        WITH last(nodes(path)) AS b
        MATCH (b)<-[e:FRAG]-(c:Mol)
        WHERE e.prop_pharmfp IS NOT NULL
        WITH %(metric)s(e.prop_pharmfp, '%(vector)s') AS sim, c.smiles AS smi, e.prop_synthon AS syn, c.cmpd_ids AS ids
        WHERE sim >= %(threshold)d
        AND NOT e.prop_synthon='%(synthon)s'
        RETURN smi, syn, sim, ids
        """ % {
            "smiles": smiles,
            "synthon": synthon,
            "num_hops": num_hops,
            "metric": metric,
            "vector": vector,
            "threshold": threshold
        }
    else:
        query = """
        MATCH (a:F2 {smiles: '%(smiles)s'})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
        WHERE e.prop_pharmfp IS NOT NULL
        WITH %(metric)s(e.prop_pharmfp, '%(vector)s') AS sim, c.smiles AS smi, e.prop_synthon AS syn, c.cmpd_ids AS ids
        WHERE sim >= %(threshold)d
        AND NOT e.prop_synthon='%(synthon)s'
        RETURN smi, syn, sim, ids
        """ % {
            "smiles": smiles,
            "synthon": synthon,
            "num_hops": num_hops,
            "metric": metric,
            "vector": vector,
            "threshold": threshold
        }

    if limit:
        query = query + f" LIMIT {limit}"

    try:
        records = run_query(query, timeout=timeout)
    except Exception as e:
        mrich.error(index, e)
        raise Exception(f"{smiles=} {synthon=} {e}")

    results = []
    for record in records:
        results.append(
            (
                record["smi"],  # expansion
                record["syn"],  # synthon
                record["sim"],  # similarity
                record["ids"],  # compound names / IDs
            )
        )

    if cache_dir:
        json.dump(results, open(cache_file, "wt"), indent=2)

    logging.info(f"Success {index} {smiles} {synthon} #results: {len(results)}")

    return results

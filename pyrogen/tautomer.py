#!/usr/bin/env python
# -*- coding: utf-8 -*-
import copy
from itertools import tee
import logging
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, BondStereo, BondDir


__author__ = 'Matt Swain'
__email__ = 'm.swain@me.com'
__license__ = 'MIT'

log = logging.getLogger('tautomer')

BONDMAP = {'-': BondType.SINGLE, '=': BondType.DOUBLE, '#': BondType.TRIPLE, ':': BondType.AROMATIC}
CHARGEMAP = {'+': 1, '0': 0, '-': -1}

tautomer_transforms = [
    {'name': '1,3 (thio)keto/enol f', 'smarts': '[CX4!H0][C]=[O,S,Se,Te;X1]'},
    {'name': '1,3 (thio)keto/enol r', 'smarts': '[O,S,Se,Te;X2!H0][C]=[C]'},
    {'name': '1,5 (thio)keto/enol f', 'smarts': '[CX4,NX3;!H0][C]=[C][CH0]=[O,S,Se,Te;X1]'},
    {'name': '1,5 (thio)keto/enol r', 'smarts': '[O,S,Se,Te;X2!H0][CH0]=,:[C][C]=,:[C,N]'},
    {'name': 'aliphatic imine f', 'smarts': '[CX4!H0][C]=[NX2]'},
    {'name': 'aliphatic imine r', 'smarts': '[NX3!H0][C]=[CX3]'},
    {'name': 'special imine f', 'smarts': '[N!H0][C]=[CX3R0]'},
    {'name': 'special imine r', 'smarts': '[CX4!H0][c]=,:[n]'},
    {'name': '1,3 aromatic heteroatom H shift f', 'smarts': '[#7!H0][#6R1]=[O,#7X2]'},
    {'name': '1,3 aromatic heteroatom H shift r', 'smarts': '[O,#7;!H0][#6R1]=,:[#7X2]'},
    {'name': '1,3 heteroatom H shift', 'smarts': '[#7,S,O,Se,Te;!H0][#7X2,#6,#15]=[#7,#16,#8,Se,Te]'},
    {'name': '1,5 aromatic heteroatom H shift', 'smarts': '[n,s,o;!H0]:[c,n]:[c]:[c,n]:[n,s,o;H0]'},
    {'name': '1,5 aromatic heteroatom H shift f', 'smarts': '[#7,#16,#8,Se,Te;!H0][#6,nX2]=,:[#6,nX2][#6,#7X2]=,:[#7X2,S,O,Se,Te]'},
    {'name': '1,5 aromatic heteroatom H shift r', 'smarts': '[#7,S,O,Se,Te;!H0][#6,#7X2]=,:[#6,nX2][#6,nX2]=,:[#7,#16,#8,Se,Te]'},
    {'name': '1,7 aromatic heteroatom H shift f', 'smarts': '[#7,#8,#16,Se,Te;!H0][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6][#6,#7X2]=,:[#7X2,S,O,Se,Te,CX3]'},
    {'name': '1,7 aromatic heteroatom H shift r', 'smarts': '[#7,S,O,Se,Te,CX4;!H0][#6,#7X2]=,:[#6][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[NX2,S,O,Se,Te]'},
    {'name': '1,9 aromatic heteroatom H shift f', 'smarts': '[#7,O;!H0][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#7,O]'},
    {'name': '1,11 aromatic heteroatom H shift f', 'smarts': '[#7,O;!H0][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#7X2,O]'},
    {'name': 'furanone f', 'smarts': '[O,S,N;!H0][#6X3r5;$([#6][!#6])]=,:[#6X3r5]'},
    {'name': 'furanone r', 'smarts': '[#6r5!H0][#6X3r5;$([#6][!#6])]=[O,S,N]'},
    {'name': 'keten/ynol f', 'smarts': '[C!H0]=[C]=[O,S,Se,Te;X1]', 'bonds': '#-'},
    {'name': 'keten/ynol r', 'smarts': '[O,S,Se,Te;!H0X2][C]#[C]', 'bonds': '=='},
    {'name': 'ionic nitro/aci-nitro f', 'smarts': '[C!H0][N+;$([N][O-])]=[O]'},
    {'name': 'ionic nitro/aci-nitro r', 'smarts': '[O!H0][N+;$([N][O-])]=[C]'},
    {'name': 'oxim/nitroso f', 'smarts': '[O!H0][N]=[C]'},
    {'name': 'oxim/nitroso r', 'smarts': '[C!H0][N]=[O]'},
    {'name': 'oxim/nitroso via phenol f', 'smarts': '[O!H0][N]=[C][C]=[C][C]=[OH0]'},
    {'name': 'oxim/nitroso via phenol r', 'smarts': '[O!H0][c]:[c]:[c]:[c][N]=[OH0]'},
    {'name': 'cyano/iso-cyanic acid f', 'smarts': '[O!H0][C]#[N]', 'bonds': '=='},
    {'name': 'cyano/iso-cyanic acid r', 'smarts': '[N!H0]=[C]=[O]', 'bonds': '#-'},
    {'name': 'formamidinesulfinic acid f', 'smarts': '[O,N;!H0][C][S,Se,Te]=[O]', 'bonds': '=--'},
    {'name': 'formamidinesulfinic acid r', 'smarts': '[O!H0][S,Se,Te][C]=[O,N]', 'bonds': '=--'},
    {'name': 'isocyanide f', 'smarts': '[C-0!H0]#[N+0]', 'bonds': '#', 'charges': '-+'},
    {'name': 'isocyanide r', 'smarts': '[N+!H0]#[C-]', 'bonds': '#', 'charges': '-+'},
    {'name': 'phosphonic acid f', 'smarts': '[OH][PH0]', 'bonds': '='},
    {'name': 'phosphonic acid r', 'smarts': '[PH]=[O]', 'bonds': '-'}
]
for transform in tautomer_transforms:
    transform['smarts'] = Chem.MolFromSmarts(transform['smarts'].encode('utf8'))

tautomer_scores = [
    {'name': 'benzoquinone', 'smarts': '[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]', 'score': 25},
    {'name': 'oxim', 'smarts': '[#6]=[N][OH]', 'score': 4},
    {'name': 'C=O', 'smarts': '[#6]=,:[#8]', 'score': 2},
    {'name': 'N=O', 'smarts': '[#7]=,:[#8]', 'score': 2},
    {'name': 'P=O', 'smarts': '[#15]=,:[#8]', 'score': 2},
    {'name': 'C=hetero', 'smarts': '[#6]=[!#1;!#6]', 'score': 1},
    {'name': 'methyl', 'smarts': '[CX4H3]', 'score': 1},
    {'name': 'guanidine terminal=N', 'smarts': '[#7][#6](=[NR0])[#7H0]', 'score': 1},
    {'name': 'guanidine endocyclic=N', 'smarts': '[#7;R][#6;R]([N])=[#7;R]', 'score': 2},
    {'name': 'aci-nitro', 'smarts': '[#6]=[N+]([O-])[OH]', 'score': -4},
]
for tscore in tautomer_scores:
    tscore['smarts'] = Chem.MolFromSmarts(tscore['smarts'])

def tautomer_score(mol):

    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    log.debug('Tautomer: %s', smiles)
    score = 0
    # Add aromatic ring scores
    ssr = Chem.GetSymmSSSR(mol)
    for ring in ssr:
        btypes = {mol.GetBondBetweenAtoms(*pair).GetBondType() for pair in _pairwise(ring)}
        elements = {mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring}
        if btypes == {BondType.AROMATIC}:
            log.debug('Score +100 (aromatic ring)')
            score += 100
            if elements == {6}:
                log.debug('Score +150 (carbocyclic aromatic ring)')
                score += 150
    # Add SMARTS scores
    for tscore in tautomer_scores:
        for match in mol.GetSubstructMatches(tscore['smarts']):
            log.debug('Score %+d (%s)', tscore['score'], tscore['name'])
            score += tscore['score']
    # Add (P,S,Se,Te)-H scores
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in {15, 16, 34, 52}:
            hs = atom.GetTotalNumHs()
            if hs:
                log.debug('Score %+d (%s-H bonds)', -hs, atom.GetSymbol())
                score -= hs

    return score


def canonical_tautomer(mol, max_tautomers=1000):
    """Enumerate all possible tautomers and return a canonical tautomer based on a scoring system.

    :param mol: An RDKit Mol object.
    :param max_tautomers: The maximum number of tautomers to enumerate (limit to prevent combinatorial explosion)
    """
    tautomers = enumerate_tautomers(mol, max_tautomers)
    if len(tautomers) == 1:
        return tautomers[0]
    # Calculate score for each tautomer
    highest = None
    for t in tautomers:
        smiles = Chem.MolToSmiles(t, isomericSmiles=True)
        log.debug('Tautomer: %s', smiles)
        score = 0
        # Add aromatic ring scores
        ssr = Chem.GetSymmSSSR(t)
        for ring in ssr:
            btypes = {t.GetBondBetweenAtoms(*pair).GetBondType() for pair in _pairwise(ring)}
            elements = {t.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring}
            if btypes == {BondType.AROMATIC}:
                log.debug('Score +100 (aromatic ring)')
                score += 100
                if elements == {6}:
                    log.debug('Score +150 (carbocyclic aromatic ring)')
                    score += 150
        # Add SMARTS scores
        for tscore in tautomer_scores:
            for match in t.GetSubstructMatches(tscore['smarts']):
                log.debug('Score %+d (%s)', tscore['score'], tscore['name'])
                score += tscore['score']
        # Add (P,S,Se,Te)-H scores
        for atom in t.GetAtoms():
            if atom.GetAtomicNum() in {15, 16, 34, 52}:
                hs = atom.GetTotalNumHs()
                if hs:
                    log.debug('Score %+d (%s-H bonds)', -hs, atom.GetSymbol())
                    score -= hs
        # Set as highest if score higher or if score equal and smiles comes first alphabetically
        if not highest or highest['score'] < score or (highest['score'] == score and smiles < highest['smiles']):
            log.debug('New highest tautomer: %s (%s)', smiles, score)
            highest = {'smiles': smiles, 'tautomer': t, 'score': score}
    return highest['tautomer']


def enumerate_tautomers(mol, max_tautomers=1000):
    """Enumerate all possible tautomers and return them as a list.

    :param mol: An RDKit Mol object.
    :param max_tautomers: The maximum number of tautomers to enumerate (limit to prevent combinatorial explosion)
    """
    tautomers = {Chem.MolToSmiles(mol, isomericSmiles=True): copy.deepcopy(mol)}
    done = set()
    while len(tautomers) < max_tautomers:
        for tsmiles in sorted(tautomers):
            if tsmiles in done:
                continue
            for transform in tautomer_transforms:
                for match in tautomers[tsmiles].GetSubstructMatches(transform['smarts']):
                    # Adjust hydrogens
                    product = copy.deepcopy(tautomers[tsmiles])
                    first = product.GetAtomWithIdx(match[0])
                    last = product.GetAtomWithIdx(match[-1])
                    first.SetNumExplicitHs(max(0, first.GetNumExplicitHs() - 1))
                    last.SetNumExplicitHs(last.GetTotalNumHs() + 1)
                    # Adjust bond orders
                    for bi, pair in enumerate(_pairwise(match)):
                        if 'bonds' in transform:
                            product.GetBondBetweenAtoms(*pair).SetBondType(BONDMAP[transform['bonds'][bi]])
                        else:
                            product.GetBondBetweenAtoms(*pair).SetBondType(BondType.DOUBLE if bi % 2 == 0 else BondType.SINGLE)
                    # Adjust charges
                    if 'charges' in transform:
                        for ci, idx in enumerate(match):
                            atom = product.GetAtomWithIdx(idx)
                            atom.SetFormalCharge(atom.GetFormalCharge() + CHARGEMAP[transform['charges'][ci]])
                    try:
                        Chem.SanitizeMol(product)
                        smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                        log.debug('Applied rule: %s to %s', transform['name'], tsmiles)
                        if smiles not in tautomers:
                            log.debug('New tautomer produced: %s' % smiles)
                            tautomers[smiles] = product
                        else:
                            log.debug('Previous tautomer produced again: %s' % smiles)
                    except ValueError:
                        log.debug('ValueError')
            done.add(tsmiles)
        if len(tautomers) == len(done):
            break
    else:
        log.warn('Tautomer enumeration stopped at maximum %s', max_tautomers)
    # Clean up stereochemistry
    for tautomer in tautomers.values():
        Chem.AssignStereochemistry(tautomer, force=True, cleanIt=True)
        for bond in tautomer.GetBonds():
            if bond.GetBondType() == BondType.DOUBLE and bond.GetStereo() > BondStereo.STEREOANY:
                begin = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                for othertautomer in tautomers.values():
                    if not othertautomer.GetBondBetweenAtoms(begin, end).GetBondType() == BondType.DOUBLE:
                        neighbours = tautomer.GetAtomWithIdx(begin).GetBonds() + tautomer.GetAtomWithIdx(end).GetBonds()
                        for otherbond in neighbours:
                            if otherbond.GetBondDir() in {BondDir.ENDUPRIGHT, BondDir.ENDDOWNRIGHT}:
                                otherbond.SetBondDir(BondDir.NONE)
                        Chem.AssignStereochemistry(tautomer, force=True, cleanIt=True)
                        log.debug('Removed stereochemistry from unfixed double bond')
                        break
    return tautomers.values()


def _pairwise(iterable):
    """Utility function to iterate in a pairwise fashion."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


tautomer_enumeration_tests = [
    ('1,3 keto/enol tautomer', 'C1(=CCCCC1)O', {'OC1=CCCCC1', 'O=C1CCCCC1'}),
    ('1,3 keto/enol tautomer', 'C1(CCCCC1)=O', {'OC1=CCCCC1', 'O=C1CCCCC1'}),
    ('Acetophenone keto/enol tautomer', 'C(=C)(O)C1=CC=CC=C1', {'C=C(O)c1ccccc1', 'CC(=O)c1ccccc1'}),
    ('Acetone keto/enol tautomer', 'CC(C)=O', {'CC(C)=O', 'C=C(C)O'}),
    ('keto/enol tautomer', 'OC(C)=C(C)C', {'C=C(O)C(C)C', 'CC(C)=C(C)O', 'CC(=O)C(C)C'}),
    ('1-phenyl-2-propanone enol/keto', 'c1(ccccc1)CC(=O)C', {'C=C(O)Cc1ccccc1', 'CC(=O)Cc1ccccc1', 'CC(O)=Cc1ccccc1'}),
    ('1,5 keto/enol tautomer', 'Oc1nccc2cc[nH]c(=N)c12', {'Nc1nccc2ccnc(O)c12', 'N=c1nccc2cc[nH]c(O)c1-2', 'Nc1[nH]ccc2ccnc(=O)c1-2', 'N=c1[nH]ccc2ccnc(O)c21', 'Nc1nccc2cc[nH]c(=O)c12', 'N=c1[nH]ccc2cc[nH]c(=O)c21'}),
    ('1,5 keto/enol tautomer', 'C1(C=CCCC1)=O', {'O=C1C=CCCC1', 'OC1=CCC=CC1', 'OC1=CC=CCC1', 'O=C1CC=CCC1', 'OC1=CCCC=C1'}),
    ('1,5 keto/enol tautomer', 'C1(=CC=CCC1)O', {'O=C1C=CCCC1', 'OC1=CCC=CC1', 'OC1=CC=CCC1', 'O=C1CC=CCC1', 'OC1=CCCC=C1'}),
    ('aliphatic imine tautomer', 'C1(CCCCC1)=N', {'N=C1CCCCC1', 'NC1=CCCCC1'}),
    ('aliphatic imine tautomer', 'C1(=CCCCC1)N', {'N=C1CCCCC1', 'NC1=CCCCC1'}),
    ('special imine tautomer', 'C1(C=CC=CN1)=CC', {'CC=C1C=CC=CN1', 'CCc1ccccn1', 'CC=C1C=CCC=N1'}),
    ('special imine tautomer', 'C1(=NC=CC=C1)CC', {'CC=C1C=CC=CN1', 'CCc1ccccn1', 'CC=C1C=CCC=N1'}),
    ('1,3 aromatic heteroatom H shift', 'O=c1cccc[nH]1', {'Oc1ccccn1', 'O=c1cccc[nH]1'}),
    ('1,3 aromatic heteroatom H shift', 'Oc1ccccn1', {'Oc1ccccn1', 'O=c1cccc[nH]1'}),
    ('1,3 aromatic heteroatom H shift', 'Oc1ncc[nH]1', {'Oc1ncc[nH]1', 'O=c1[nH]cc[nH]1'}),
    ('1,3 heteroatom H shift', 'OC(C)=NC', {'CN=C(C)O', 'CNC(C)=O', 'C=C(O)NC'}),
    ('1,3 heteroatom H shift', 'CNC(C)=O', {'CN=C(C)O', 'CNC(C)=O', 'C=C(O)NC'}),
    ('1,3 heteroatom H shift', 'S=C(N)N', {'N=C(N)S', 'NC(N)=S'}),
    ('1,3 heteroatom H shift', 'SC(N)=N', {'N=C(N)S', 'NC(N)=S'}),
    ('1,3 heteroatom H shift', 'N=c1[nH]ccn(C)1', {'Cn1ccnc1N', 'Cn1cc[nH]c1=N'}),
    ('1,3 heteroatom H shift', 'CN=c1[nH]cncc1', {'CN=c1ccnc[nH]1', 'CNc1ccncn1', 'CN=c1cc[nH]cn1'}),
    ('1,5 aromatic heteroatom H shift', 'Oc1cccc2ccncc12', {'O=c1cccc2cc[nH]cc1-2', 'Oc1cccc2ccncc12'}),
    ('1,5 aromatic heteroatom H shift', 'O=c1cccc2cc[nH]cc1-2', {'O=c1cccc2cc[nH]cc1-2', 'Oc1cccc2ccncc12'}),
    ('1,5 aromatic heteroatom H shift', 'Cc1n[nH]c2ncnn12', {'C=C1NNc2ncnn21', 'Cc1n[nH]c2ncnn12', 'Cc1nnc2[nH]cnn12', 'C=C1NN=C2N=CNN12', 'Cc1nnc2nc[nH]n12', 'C=C1NN=C2NC=NN12'}),
    ('1,5 aromatic heteroatom H shift', 'Cc1nnc2nc[nH]n12', {'C=C1NNc2ncnn21', 'Cc1n[nH]c2ncnn12', 'Cc1nnc2[nH]cnn12', 'C=C1NN=C2N=CNN12', 'Cc1nnc2nc[nH]n12', 'C=C1NN=C2NC=NN12'}),
    ('1,5 aromatic heteroatom H shift', 'Oc1ccncc1', {'Oc1ccncc1', 'O=c1cc[nH]cc1'}),
    ('1,5 aromatic heteroatom H shift', 'Oc1c(cccc3)c3nc2ccncc12', {'Oc1c2ccccc2nc2ccncc12', 'O=c1c2ccccc2nc2cc[nH]cc1-2', 'O=c1c2ccccc2[nH]c2ccncc21'}),
    ('1,3 and 1,5 aromatic heteroatom H shift', 'Oc1ncncc1', {'Oc1ccncn1', 'O=c1ccnc[nH]1', 'O=c1cc[nH]cn1'}),
    ('1,5 aromatic heteroatom H shift', 'C2(=C1C(=NC=N1)[NH]C(=N2)N)O', {'N=c1[nH]c2ncnc-2c(O)[nH]1', 'Nc1nc2nc[nH]c2c(O)n1', 'N=c1nc(O)c2nc[nH]c2[nH]1', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1', 'Nc1nc2ncnc-2c(O)[nH]1', 'N=c1nc2nc[nH]c2c(O)[nH]1', 'N=c1nc(O)c2[nH]cnc2[nH]1', 'Nc1nc(O)c2ncnc-2[nH]1', 'Nc1nc(=O)c2nc[nH]c2[nH]1', 'Nc1nc(=O)c2[nH]cnc2[nH]1', 'Nc1nc2[nH]cnc2c(O)n1', 'N=c1nc2[nH]cnc2c(O)[nH]1', 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'Nc1nc2nc[nH]c2c(=O)[nH]1', 'N=c1[nH]c2nc[nH]c2c(=O)[nH]1'}),
    ('1,5 aromatic heteroatom H shift', 'C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O', {'N=c1[nH]c2ncnc-2c(O)[nH]1', 'Nc1nc2nc[nH]c2c(O)n1', 'N=c1nc(O)c2nc[nH]c2[nH]1', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1', 'Nc1nc2ncnc-2c(O)[nH]1', 'N=c1nc2nc[nH]c2c(O)[nH]1', 'N=c1nc(O)c2[nH]cnc2[nH]1', 'Nc1nc(O)c2ncnc-2[nH]1', 'Nc1nc(=O)c2nc[nH]c2[nH]1', 'Nc1nc(=O)c2[nH]cnc2[nH]1', 'Nc1nc2[nH]cnc2c(O)n1', 'N=c1nc2[nH]cnc2c(O)[nH]1', 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'Nc1nc2nc[nH]c2c(=O)[nH]1', 'N=c1[nH]c2nc[nH]c2c(=O)[nH]1'}),
    ('1,5 aromatic heteroatom H shift', 'Oc1n(C)ncc1', {'Cn1nccc1O', 'CN1N=CCC1=O', 'Cn1[nH]ccc1=O'}),
    ('1,5 aromatic heteroatom H shift', 'O=c1nc2[nH]ccn2cc1', {'O=c1ccn2cc[nH]c2n1', 'Oc1ccn2ccnc2n1', 'O=c1ccn2ccnc2[nH]1'}),
    ('1,5 aromatic heteroatom H shift', 'N=c1nc[nH]cc1', {'N=c1cc[nH]cn1', 'N=c1ccnc[nH]1', 'Nc1ccncn1'}),
    ('1,5 aromatic heteroatom H shift', 'N=c(c1)ccn2cc[nH]c12', {'N=c1ccn2cc[nH]c2c1', 'Nc1ccn2ccnc2c1'}),
    ('1,5 aromatic heteroatom H shift', 'CN=c1nc[nH]cc1', {'CN=c1ccnc[nH]1', 'CNc1ccncn1', 'CN=c1cc[nH]cn1'}),
    ('1,7 aromatic heteroatom H shift', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', {'c1ccc2c(c1)=NC(C1=NC3C=CC=CC3=N1)N=2', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', 'c1ccc2[nH]c(C3=NC4C=CC=CC4=N3)nc2c1', 'c1ccc2[nH]c(C3N=c4ccccc4=N3)nc2c1', 'c1ccc2c(c1)=NC(=C1N=C3C=CC=CC3N1)N=2', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'}),
    ('1,7 aromatic heteroatom H shift', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2', {'c1ccc2c(c1)=NC(C1=NC3C=CC=CC3=N1)N=2', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', 'c1ccc2[nH]c(C3=NC4C=CC=CC4=N3)nc2c1', 'c1ccc2[nH]c(C3N=c4ccccc4=N3)nc2c1', 'c1ccc2c(c1)=NC(=C1N=C3C=CC=CC3N1)N=2', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'}),
    ('1,9 aromatic heteroatom H shift', 'CNc1ccnc2ncnn21', {'CN=c1cc[nH]c2ncnn21', 'CN=c1ccnc2[nH]cnn21', 'CN=c1ccnc2nc[nH]n21', 'CNc1ccnc2ncnn21'}),
    ('1,9 aromatic heteroatom H shift', 'CN=c1ccnc2nc[nH]n21', {'CN=c1cc[nH]c2ncnn21', 'CN=c1ccnc2[nH]cnn21', 'CN=c1ccnc2nc[nH]n21', 'CNc1ccnc2ncnn21'}),
    ('1,11 aromatic heteroatom H shift', 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', {'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', 'N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', 'N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1'}),
    ('1,11 aromatic heteroatom H shift', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', {'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', 'N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', 'N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1'}),
    ('heterocyclic tautomer', 'n1ccc2ccc[nH]c12', {'c1cc2cccnc2[nH]1', 'c1cc2ccc[nH]c-2n1'}),
    ('heterocyclic tautomer', 'c1cc(=O)[nH]c2nccn12', {'O=c1ccn2cc[nH]c2n1', 'Oc1ccn2ccnc2n1', 'O=c1ccn2ccnc2[nH]1'}),
    ('heterocyclic tautomer', 'c1cnc2c[nH]ccc12', {'c1cc2cc[nH]cc-2n1', 'c1cc2ccncc2[nH]1'}),
    ('heterocyclic tautomer', 'n1ccc2c[nH]ccc12', {'c1cc2cnccc2[nH]1', 'c1cc2c[nH]ccc-2n1'}),
    ('heterocyclic tautomer', 'c1cnc2ccc[nH]c12', {'c1cc2[nH]cccc-2n1', 'c1cc2ncccc2[nH]1'}),
    ('furanone tautomer', 'C1=CC=C(O1)O', {'Oc1ccco1', 'O=C1CC=CO1'}),
    ('furanone tautomer', 'O=C1CC=CO1', {'Oc1ccco1', 'O=C1CC=CO1'}),
    ('keten/ynol tautomer', 'CC=C=O', {'CC=C=O', 'CC#CO'}),
    ('keten/ynol tautomer', 'CC#CO', {'CC=C=O', 'CC#CO'}),
    ('ionic nitro/aci-nitro tautomer', 'C([N+](=O)[O-])C', {'CC[N+](=O)[O-]', 'CC=[N+]([O-])O'}),
    ('ionic nitro/aci-nitro tautomer', 'C(=[N+](O)[O-])C', {'CC[N+](=O)[O-]', 'CC=[N+]([O-])O'}),
    ('oxim nitroso tautomer', 'CC(C)=NO', {'CC(C)N=O', 'CC(C)=NO', 'C=C(C)NO'}),
    ('oxim nitroso tautomer', 'CC(C)N=O', {'CC(C)N=O', 'CC(C)=NO', 'C=C(C)NO'}),
    ('oxim/nitroso tautomer via phenol', 'O=Nc1ccc(O)cc1', {'O=NC1C=CC(=O)C=C1', 'O=C1C=CC(=NO)C=C1', 'O=Nc1ccc(O)cc1'}),
    ('oxim/nitroso tautomer via phenol', 'O=C1C=CC(=NO)C=C1', {'O=NC1C=CC(=O)C=C1', 'O=C1C=CC(=NO)C=C1', 'O=Nc1ccc(O)cc1'}),
    ('cyano/iso-cyanic acid tautomer', 'C(#N)O', {'N#CO', 'N=C=O'}),
    ('cyano/iso-cyanic acid tautomer', 'C(=N)=O', {'N#CO', 'N=C=O'}),
    ('isocyanide tautomer', 'C#N', {'[C-]#[NH+]', 'C#N'}),
    ('isocyanide tautomer', '[C-]#[NH+]', {'[C-]#[NH+]', 'C#N'}),
    ('Remove stereochemistry from mobile double bonds', 'c1(ccccc1)/C=C(/O)\\C', {'C=C(O)Cc1ccccc1', 'CC(O)=Cc1ccccc1', 'CC(=O)Cc1ccccc1'}),
    ('Remove stereochemistry from mobile double bonds', 'C/C=C/C(C)=O', {'C=C(O)C=CC', 'C=CCC(=C)O', 'CC=CC(C)=O', 'C=CCC(C)=O', 'C=CC=C(C)O'}),
    ('Remove stereochemistry from mobile double bonds', 'C/C=C\C(C)=O', {'C=C(O)C=CC', 'C=CCC(=C)O', 'CC=CC(C)=O', 'C=CCC(C)=O', 'C=CC=C(C)O'}),
]

tautomer_canonicalization_tests = [
    ('1,3 keto/enol tautomer', 'C1(=CCCCC1)O', 'O=C1CCCCC1'),
    ('1,3 keto/enol tautomer', 'C1(CCCCC1)=O', 'O=C1CCCCC1'),
    ('Acetophenone keto/enol tautomer', 'C(=C)(O)C1=CC=CC=C1', 'CC(=O)c1ccccc1'),
    ('Acetone keto/enol tautomer', 'CC(C)=O', 'CC(C)=O'),
    ('keto/enol tautomer', 'OC(C)=C(C)C', 'CC(=O)C(C)C'),
    ('1-phenyl-2-propanone enol/keto', 'c1(ccccc1)CC(=O)C', 'CC(=O)Cc1ccccc1'),
    ('1,5 keto/enol tautomer', 'Oc1nccc2cc[nH]c(=N)c12', 'N=c1[nH]ccc2cc[nH]c(=O)c21'),
    ('1,5 keto/enol tautomer', 'C1(C=CCCC1)=O', 'O=C1C=CCCC1'),
    ('1,5 keto/enol tautomer', 'C1(=CC=CCC1)O', 'O=C1C=CCCC1'),
    ('aliphatic imine tautomer', 'C1(CCCCC1)=N', 'N=C1CCCCC1'),
    ('aliphatic imine tautomer', 'C1(=CCCCC1)N', 'N=C1CCCCC1'),
    ('special imine tautomer', 'C1(C=CC=CN1)=CC', 'CCc1ccccn1'),
    ('special imine tautomer', 'C1(=NC=CC=C1)CC', 'CCc1ccccn1'),
    ('1,3 aromatic heteroatom H shift', 'O=c1cccc[nH]1', 'O=c1cccc[nH]1'),
    ('1,3 aromatic heteroatom H shift', 'Oc1ccccn1', 'O=c1cccc[nH]1'),
    ('1,3 aromatic heteroatom H shift', 'Oc1ncc[nH]1', 'O=c1[nH]cc[nH]1'),
    ('1,3 heteroatom H shift', 'OC(C)=NC', 'CNC(C)=O'),
    ('1,3 heteroatom H shift', 'CNC(C)=O', 'CNC(C)=O'),
    ('1,3 heteroatom H shift', 'S=C(N)N', 'NC(N)=S'),
    ('1,3 heteroatom H shift', 'SC(N)=N', 'NC(N)=S'),
    ('1,3 heteroatom H shift', 'N=c1[nH]ccn(C)1', 'Cn1cc[nH]c1=N'),
    ('1,3 heteroatom H shift', 'CN=c1[nH]cncc1', 'CN=c1cc[nH]cn1'),
    ('1,5 aromatic heteroatom H shift', 'Oc1cccc2ccncc12', 'Oc1cccc2ccncc12'),
    ('1,5 aromatic heteroatom H shift', 'O=c1cccc2cc[nH]cc1-2', 'Oc1cccc2ccncc12'),
    ('1,5 aromatic heteroatom H shift', 'Cc1n[nH]c2ncnn12', 'Cc1n[nH]c2ncnn12'),
    ('1,5 aromatic heteroatom H shift', 'Cc1nnc2nc[nH]n12', 'Cc1n[nH]c2ncnn12'),
    ('1,5 aromatic heteroatom H shift', 'Oc1ccncc1', 'O=c1cc[nH]cc1'),
    ('1,5 aromatic heteroatom H shift', 'Oc1c(cccc3)c3nc2ccncc12', 'O=c1c2ccccc2[nH]c2ccncc21'),
    ('1,3 and 1,5 aromatic heteroatom H shift', 'Oc1ncncc1', 'O=c1cc[nH]cn1'),
    ('1,5 aromatic heteroatom H shift', 'C2(=C1C(=NC=N1)[NH]C(=N2)N)O', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1'),
    ('1,5 aromatic heteroatom H shift', 'C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1'),
    ('1,5 aromatic heteroatom H shift', 'Oc1n(C)ncc1', 'Cn1[nH]ccc1=O'),
    ('1,5 aromatic heteroatom H shift', 'O=c1nc2[nH]ccn2cc1', 'O=c1ccn2cc[nH]c2n1'),
    ('1,5 aromatic heteroatom H shift', 'N=c1nc[nH]cc1', 'N=c1cc[nH]cn1'),
    ('1,5 aromatic heteroatom H shift', 'N=c(c1)ccn2cc[nH]c12', 'N=c1ccn2cc[nH]c2c1'),
    ('1,5 aromatic heteroatom H shift', 'CN=c1nc[nH]cc1', 'CN=c1cc[nH]cn1'),
    ('1,7 aromatic heteroatom H shift', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1'),
    ('1,7 aromatic heteroatom H shift', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1'),
    ('1,9 aromatic heteroatom H shift', 'CNc1ccnc2ncnn21', 'CN=c1cc[nH]c2ncnn21'),
    ('1,9 aromatic heteroatom H shift', 'CN=c1ccnc2nc[nH]n21', 'CN=c1cc[nH]c2ncnn21'),
    ('1,11 aromatic heteroatom H shift', 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1'),
    ('1,11 aromatic heteroatom H shift', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1'),
    ('heterocyclic tautomer', 'n1ccc2ccc[nH]c12', 'c1cc2cccnc2[nH]1'),
    ('heterocyclic tautomer', 'c1cc(=O)[nH]c2nccn12', 'O=c1ccn2cc[nH]c2n1'),
    ('heterocyclic tautomer', 'c1cnc2c[nH]ccc12', 'c1cc2ccncc2[nH]1'),
    ('heterocyclic tautomer', 'n1ccc2c[nH]ccc12', 'c1cc2cnccc2[nH]1'),
    ('heterocyclic tautomer', 'c1cnc2ccc[nH]c12', 'c1cc2ncccc2[nH]1'),
    ('furanone tautomer', 'C1=CC=C(O1)O', 'Oc1ccco1'),
    ('furanone tautomer', 'O=C1CC=CO1', 'Oc1ccco1'),
    ('keten/ynol tautomer', 'CC=C=O', 'CC=C=O'),
    ('keten/ynol tautomer', 'CC#CO', 'CC=C=O'),
    ('ionic nitro/aci-nitro tautomer', 'C([N+](=O)[O-])C', 'CC[N+](=O)[O-]'),
    ('ionic nitro/aci-nitro tautomer', 'C(=[N+](O)[O-])C', 'CC[N+](=O)[O-]'),
    ('oxim nitroso tautomer', 'CC(C)=NO', 'CC(C)=NO'),
    ('oxim nitroso tautomer', 'CC(C)N=O', 'CC(C)=NO'),
    ('oxim/nitroso tautomer via phenol', 'O=Nc1ccc(O)cc1', 'O=Nc1ccc(O)cc1'),
    ('oxim/nitroso tautomer via phenol', 'O=C1C=CC(=NO)C=C1', 'O=Nc1ccc(O)cc1'),
    ('cyano/iso-cyanic acid tautomer', 'C(#N)O', 'N=C=O'),
    ('cyano/iso-cyanic acid tautomer', 'C(=N)=O', 'N=C=O'),
    ('formamidinesulfinic acid tautomer', '[S](=O)(=O)C(N)N', 'N=C(N)S(=O)O'),
    ('formamidinesulfinic acid tautomer', '[S](=O)(O)C(=N)N', 'N=C(N)S(=O)O'),
    ('isocyanide tautomer', 'C#N', 'C#N'),
    ('isocyanide tautomer', '[C-]#[NH+]', 'C#N'),
]

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    for desc, smiles, tautomers in tautomer_enumeration_tests:
        mol = Chem.MolFromSmiles(smiles)
        Chem.SanitizeMol(mol)
        assert {Chem.MolToSmiles(t, isomericSmiles=True) for t in enumerate_tautomers(mol)} == tautomers
        log.info('%s => %s', smiles, tautomers)

    for desc, smiles, tautomer in tautomer_canonicalization_tests:
        mol = Chem.MolFromSmiles(smiles)
        Chem.SanitizeMol(mol)
        assert Chem.MolToSmiles(canonical_tautomer(mol), isomericSmiles=True) == tautomer
        log.info('%s => %s', smiles, tautomer)

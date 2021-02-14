import mbuild as mb
import antefoyer
import numpy as np
import foyer
import parmed as pmd
import copy
import argparse
from atools.lib.chains import Alkylsilane
from atools.recipes import DualSurface, SilicaInterface, SurfaceMonolayer
from atools.lib.chains.alkylsilane_internal import Alkylsilane as AlkylsilaneInternal
from mbuild.lib.atoms import H
from foyer import Forcefield
from foyer import forcefield
from foyer.atomtyper import find_atomtypes
import simtk.openmm.app.element as elem
import foyer.element as custom_elem
import simtk.unit as u
from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app.forcefield import (NoCutoff, CutoffNonPeriodic, HBonds, AllBonds, HAngles, NonbondedGenerator, _convertParameterToNumber)


def main(args):
    length = 2
    backbone = args.backbone
    locations = [0]
    terminal_group = args.terminal_group

    cpa = AlkylsilaneInternal(chain_length=length, internal_group=backbone, locations=locations, terminal_group=terminal_group)
    h = H()
    compound = mb.Compound()
    compound.add(cpa, 'cpa')
    compound.add(h, 'h')
    mb.force_overlap(compound['cpa'], compound['cpa'].all_ports()[0], compound['h'].all_ports()[0])

    structure = cpa.to_parmed(box=None, residues=['chain'])
    ff = Forcefield(forcefield_files='../shearing-code/src/util/forcefield/oplsaa.xml')
    structure_ = ff.apply(structure)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--terminal_group', '-tg', type=str, default='methyl')
    parser.add_argument('--backbone', '-bb', type=str, default='methylene')
    args = parser.parse_args()
    
    main(args)
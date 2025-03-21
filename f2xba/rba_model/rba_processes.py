"""Implementation of RbaProcesses, RbaProcess and RbaProcessingMap classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import numpy as np
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from .rba_target_value import RbaTargetValue
from ..utils.rba_utils import get_species_refs_from_xml, get_species_refs_from_str, \
                              extract_params, translate_reaction_string


class RbaProcesses:

    def __init__(self):
        self.processes = {}
        self.processing_maps = {}

    def import_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'processes.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBAProcesses'
            self.processes = RbaProcess.import_xml(root.find('listOfProcesses'))
            self.processing_maps = RbaProcessingMap.import_xml(root.find('listOfProcessingMaps'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'processes' in m_dict:
            self.processes = RbaProcess.from_df(m_dict['processes'])
        else:
            print(f'processes not imported!')
        if 'processingMaps' in m_dict:
            self.processing_maps = RbaProcessingMap.from_df(m_dict['processingMaps'])
        else:
            print(f'processingMaps not imported!')

    def from_xba(self, rba_params, rba_model):
        """Configure Processes and Processing maps.

        Parameter functions and aggregates are created
        Cofactors are automatically added to translation map.

        'input_filter' for set 'rna'
        - comma-separated list of regular expression patterns that will be checked against the macromolecule id
            - e.g. 'mrna' to select only mrna
            - e.g. 'trna, b' to select only tRNAs and E. coli gene loci (used for rRNA)

        'input_filter' for set 'protein':
        - 'signal_peptide': select proteins having a signal peptide as per UniProt data
          plus dummy_proteins not located in the cytoplasm.
        - comma-separated list of valid RBA compartment ids (select proteins accordingly)
        - comma-separated list of regular expression patterns that will be checked against protein id

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param rba_model: rba model instance
        :type rba_model: class RbaModel
        """
        df_pmaps_data = rba_params['processing_maps']
        df_proc_data = rba_params['processes']
        df_mach_data = rba_params['machineries']
        cytoplasm_cid = rba_model.cid_mappings['cytoplasm_cid']

        # Process data
        pmap2set = {}
        for proc_id, row in df_proc_data.iterrows():
            proc_name = row['name']
            proc_set = row['set']
            proc_type = row['type']
            proc_map = row['processing_map']
            pmap2set[proc_map] = proc_set

            # add machinery, if defined
            constant = row.get('capacity_constant', np.nan)
            func = row.get('capacity_function', np.nan)
            agg = row.get('capacity_aggregate', np.nan)
            p_name = rba_model.parameters.create_parameter(f'{proc_id}_capacity', constant, func, agg)
            machinery = {}
            if p_name:
                machinery['capacity'] = RbaTargetValue.get_target_value('value', p_name)
                machinery['reactants'] = {}
                machinery['products'] = {}
                # note: df_mach_data.loc[[proc_id]] with double angular brackets to enforce dataframe being returned
                for _, mach_comp in df_mach_data.loc[[proc_id]].iterrows():
                    if mach_comp['stoic'] < 0.0:
                        machinery['reactants'][mach_comp['id']] = -mach_comp['stoic']
                    else:
                        machinery['products'][mach_comp['id']] = mach_comp['stoic']
            # determine process inputs
            inputs = []
            if proc_set == 'dna':
                inputs = list(rba_model.dna.macromolecules)

            elif proc_set == 'rna':
                if type(row['input_filter']) is not str:
                    inputs = list(rba_model.rnas.macromolecules)
                else:
                    patterns = [pattern.strip() for pattern in row['input_filter'].split(',')]
                    input_set = set()
                    for rnaid in rba_model.rnas.macromolecules:
                        for pattern in patterns:
                            if re.match(pattern, rnaid):
                                input_set.add(rnaid)
                                break
                    inputs = list(input_set)

            elif proc_set == 'protein':
                if type(row['input_filter']) is not str:
                    inputs = list(rba_model.proteins.macromolecules)
                elif row['input_filter'] == 'signal_peptide':
                    inputs = []
                    for locus in rba_model.proteins.macromolecules:
                        # add those proteins to the input, that have a signal_peptide in their UniProt data
                        if locus in rba_model.model.locus2uid:
                            p = rba_model.model.proteins[rba_model.model.locus2uid[locus]]
                            if p.has_signal_peptide is True:
                                inputs.append(locus)
                        # dummy proteins have no UniProt entries. We add them to input if they are located in cytoplasm
                        else:
                            if rba_model.proteins.macromolecules[locus].compartment != cytoplasm_cid:
                                inputs.append(locus)
                else:
                    items = {item.strip() for item in row['input_filter'].split(',')}
                    input_set = set()
                    # check if items are valid compartment ids:
                    if len(items.intersection(set(rba_model.metabolism.compartments))) == len(items):
                        for mm_id, mm in rba_model.proteins.macromolecules.items():
                            if mm.compartment in items:
                                input_set.add(mm_id)
                    else:
                        # consider items as patterns or names of macromolecules
                        for mm_id in rba_model.proteins.macromolecules:
                            for pattern in items:
                                if re.match(pattern, mm_id):
                                    input_set.add(mm_id)
                                    break
                    inputs = list(input_set)

            proc_data = {'processingMap': proc_map, 'set': proc_set, 'inputs': inputs}
            if proc_type == 'production':
                self.processes[proc_id] = RbaProcess(proc_id, name=proc_name, machinery=machinery,
                                                     productions=proc_data)
            elif proc_type == 'degradation':
                self.processes[proc_id] = RbaProcess(proc_id, name=proc_name, machinery=machinery,
                                                     degradations=proc_data)

        # processing maps
        set2components = {'dna': rba_model.dna.components, 'rna': rba_model.rnas.components,
                          'protein': rba_model.proteins.components}
        for pmapid, proc_set in pmap2set.items():
            const_proc = {}
            comp_proc = {}
            for _, row in df_pmaps_data.loc[[pmapid]].iterrows():
                mach_cost = row['machinery_cost']
                component = row['component']
                reactants, products = translate_reaction_string(row['reaction_string'])

                # constant processing data for macromolecule
                if component == 'constantProcessing':
                    const_proc['reactants'] = reactants
                    const_proc['products'] = products

                # comonent processing costs
                elif component in set2components[proc_set]:
                    # processing of a main component, e.g. a specific amino acid during translation
                    comp_proc[component] = {'cost': mach_cost, 'reactants': reactants, 'products': products}
                else:
                    # automatic creation of a set of component processings of specific type, e.g. folding of amino acids
                    for comp_id, comp in set2components[proc_set].items():
                        if comp.type == component:
                            if comp.type == 'cofactor':
                                # cofactor processing is based on cofactor consumption
                                comp_proc[comp_id] = {'cost': mach_cost, 'reactants': {comp_id: 1.0}, 'products': {}}
                            else:
                                # e.g. amino acid processing with same reaction string (e.g. folding, secretion)
                                comp_proc[comp_id] = {'cost': mach_cost, 'reactants': reactants, 'products': products}
                # configuration if individual reaction

            self.processing_maps[pmapid] = RbaProcessingMap(pmapid, const_processing=const_proc,
                                                            comp_processing=comp_proc)

        print(f'{len(self.processes):4d} processes and {len(self.processing_maps)} processing maps')

    def export_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'processes.xml')
        root = Element('RBAProcesses')

        processes = SubElement(root, 'listOfProcesses')
        for item in self.processes.values():
            processes.append(item.export_xml())

        pmaps = SubElement(root, 'listOfProcessingMaps')
        for item in self.processing_maps.values():
            pmaps.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):
        df_p = pd.DataFrame([item.to_dict() for item in self.processes.values()])
        df_p.set_index('process', inplace=True)
        data = []
        for pmid, pmap in self.processing_maps.items():
            pm_dict = pmap.to_dict()
            if ((len(pm_dict['constantProcessing'].get('reactants', {})) > 0) or
                    (len(pm_dict['constantProcessing'].get('products', {})) > 0)):
                data.append([pmid, 'constantProcessing', np.nan, pm_dict['constantProcessing']['reactants'],
                             pm_dict['constantProcessing']['products']])
            for component, comp_proc in pm_dict['componentProcessings'].items():
                data.append([pmid, component, comp_proc['machineryCost'], comp_proc['reactants'],
                             comp_proc['products']])

        df_pm = pd.DataFrame(data, columns=['processingMap', 'component', 'machineryCost',
                                            'reactants', 'products'])
        df_pm.set_index('processingMap', inplace=True)

        return {'processes': df_p, 'processingMaps': df_pm}

    def validate(self, component_ids):
        valid = True

        missing = self.ref_molecules_pmaps().difference(component_ids['species'])
        if len(missing) > 0:
            print('species used in processingMaps not defined:', missing)
            valid = False

        missing = self.ref_molecules_machinery().difference(component_ids['species']) \
            .difference(component_ids['rna']) \
            .difference(component_ids['protein'])
        if len(missing) > 0:
            print('species/macromolecules used in processes not defined:', missing)
            valid = False

        missing = self.ref_molecules_inputs().difference(component_ids['dna']) \
            .difference(component_ids['rna']) \
            .difference(component_ids['protein'])
        if len(missing) > 0:
            print('macromolecules used in processes not defined:', missing)
            valid = False

        missing = self.ref_parameters().difference(component_ids['functions']) \
            .difference(component_ids['aggregates'])
        if len(missing) > 0:
            print('function/aggregates used in aggregates not defined:', missing)
            valid = False

        return valid

    def ref_molecules_pmaps(self):
        refs = set()
        for pmap in self.processing_maps.values():
            refs |= {sid for sid in pmap.constant_processing.get('reactants', {})}
            refs |= {sid for sid in pmap.constant_processing.get('products', {})}
            for comp_proc in pmap.component_processings.values():
                refs |= {sid for sid in comp_proc['reactants']}
                refs |= {sid for sid in comp_proc['products']}
        return refs

    def ref_molecules_machinery(self):
        refs = set()
        for p in self.processes.values():
            refs |= {sid for sid in p.machinery.get('reactants', {})}
            refs |= {sid for sid in p.machinery.get('products', {})}
        return refs

    def ref_molecules_inputs(self):
        refs = set()
        for p in self.processes.values():
            refs |= {sid for sid in p.productions.get('inputs', {})}
            refs |= {sid for sid in p.degradations.get('inputs', {})}
        return refs

    def ref_molecules(self):
        refs = set()
        refs |= self.ref_molecules_pmaps()
        refs |= self.ref_molecules_machinery()
        refs |= self.ref_molecules_inputs()
        return refs

    def ref_parameters(self):
        refs = set()
        for p in self.processes.values():
            if 'capacity' in p.machinery:
                refs.add(p.machinery['capacity'].value)
        return refs


class RbaProcess:

    def __init__(self, proc_id, name=None, machinery=None, productions=None, degradations=None):
        self.id = proc_id
        self.name = name
        self.machinery = machinery if type(machinery) is dict else {}
        self.productions = productions if type(productions) is dict else {}
        self.degradations = degradations if type(degradations) is dict else {}
        self.constr_id = None

    @staticmethod
    def import_xml(processes):

        data = {}
        for process in processes.findall('process'):
            pid = process.attrib['id']
            rba_process = RbaProcess(pid)
            rba_process.name = process.get('name')
            machinery = process.find('machinery')
            if machinery is not None:
                composition = machinery.find('machineryComposition')
                reactants = get_species_refs_from_xml(composition.find('listOfReactants'))
                products = get_species_refs_from_xml(composition.find('listOfProducts'))
                capacity = RbaTargetValue.from_dict(machinery.find('capacity').attrib)
                rba_process.machinery = {'capacity': capacity, 'reactants': reactants,
                                         'products': products}
            processings = process.find('processings')
            if processings is not None:
                productions = processings.find('listOfProductions')
                if productions is not None:
                    processing = productions.find('processing')
                    processing_map = processing.attrib['processingMap']
                    set_type = processing.attrib['set']
                    lo_inputs = processing.find('listOfInputs')
                    inputs = [sref.attrib['species'] for sref in lo_inputs.findall('speciesReference')]
                    rba_process.productions = {'processingMap': processing_map, 'set': set_type, 'inputs': inputs}
                degradations = processings.find('listOfDegradations')
                if degradations is not None:
                    processing = degradations.find('processing')
                    processing_map = processing.attrib['processingMap']
                    set_type = processing.attrib['set']
                    lo_inputs = processing.find('listOfInputs')
                    inputs = [sref.attrib['species'] for sref in lo_inputs.findall('speciesReference')]
                    rba_process.degradations = {'processingMap': processing_map, 'set': set_type, 'inputs': inputs}

            data[pid] = rba_process
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for pid, row in df.iterrows():
            rba_process = RbaProcess(pid)
            rba_process.name = row['name']
            if type(row['machineryCapacity']) is str and len(row['machineryCapacity']) > 0:
                capacity_dict = extract_params(row['machineryCapacity'])
                capacity = RbaTargetValue.from_dict(capacity_dict)
                reactants = get_species_refs_from_str(row['machineryReactants'])
                products = get_species_refs_from_str(row['machineryProducts'])
                rba_process.machinery = {'capacity': capacity, 'reactants': reactants,
                                         'products': products}
            if type(row['productionProcessingMap']) is str and len(row['productionProcessingMap']) > 0:
                processing_map = row['productionProcessingMap']
                set_type = row['productionSet']
                inputs = [item.strip() for item in row['productionInputs'].split(',')]
                rba_process.productions = {'processingMap': processing_map, 'set': set_type, 'inputs': inputs}

            if type(row['degradationProcessingMap']) is str and len(row['degradationProcessingMap']) > 0:
                processing_map = row['degradationProcessingMap']
                set_type = row['degradationSet']
                inputs = [item.strip() for item in row['degradationInputs'].split(',')]
                rba_process.degradations = {'processingMap': processing_map, 'set': set_type, 'inputs': inputs}

            data[pid] = rba_process
        return data

    def export_xml(self):
        process = Element('process', {'id': self.id, 'name': self.name})

        if len(self.machinery) > 0:
            machinery = SubElement(process, 'machinery')
            composition = SubElement(machinery, 'machineryComposition')
            if len(self.machinery['reactants']) > 0:
                reactants = SubElement(composition, 'listOfReactants')
                for sid, stoic in self.machinery['reactants'].items():
                    SubElement(reactants, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
            if len(self.machinery['products']) > 0:
                products = SubElement(composition, 'listOfProducts')
                for sid, stoic in self.machinery['products'].items():
                    SubElement(products, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
            SubElement(machinery, 'capacity', self.machinery['capacity'].to_dict())
        if len(self.productions) + len(self.degradations) > 0:
            processings = SubElement(process, 'processings')
            if len(self.productions) > 0:
                productions = SubElement(processings, 'listOfProductions')
                attribs = {'processingMap': self.productions['processingMap'], 'set': self.productions['set']}
                processing = SubElement(productions, 'processing', attribs)
                inputs = SubElement(processing, 'listOfInputs')
                for sid in self.productions['inputs']:
                    SubElement(inputs, 'speciesReference', {'species': sid, 'stoichiometry': '1'})
            if len(self.degradations) > 0:
                degradations = SubElement(processings, 'listOfDegradations')
                attribs = {'processingMap': self.degradations['processingMap'], 'set': self.degradations['set']}
                processing = SubElement(degradations, 'processing', attribs)
                inputs = SubElement(processing, 'listOfInputs')
                for sid in self.degradations['inputs']:
                    SubElement(inputs, 'speciesReference', {'species': sid, 'stoichiometry': '1'})
        return process

    def to_dict(self):
        mach_capacity = ''
        mach_reactants = ''
        mach_products = ''
        if len(self.machinery) > 0:
            mach_capacity = ', '.join([f'{key}={val}'
                                       for key, val in self.machinery['capacity'].to_dict().items()])
            mach_reactants = '; '.join([f'species={species}, stoic={stoic}'
                                        for species, stoic in self.machinery['reactants'].items()])
            mach_products = '; '.join([f'species={species}, stoic={stoic}'
                                       for species, stoic in self.machinery['products'].items()])
        prod_pmap = ''
        prod_set = ''
        prod_inputs = ''
        if len(self.productions) > 0:
            prod_pmap = self.productions['processingMap']
            prod_set = self.productions['set']
            prod_inputs = ', '.join(self.productions['inputs'])

        degr_pmap = ''
        degr_set = ''
        degr_inputs = ''
        if len(self.degradations) > 0:
            degr_pmap = self.degradations['processingMap']
            degr_set = self.degradations['set']
            degr_inputs = ', '.join(self.degradations['inputs'])

        return {'process': self.id, 'name': self.name,
                'machineryCapacity': mach_capacity, 'machineryReactants': mach_reactants,
                'machineryProducts': mach_products,
                'productionProcessingMap': prod_pmap, 'productionSet': prod_set,
                'productionInputs': prod_inputs,
                'degradationProcessingMap': degr_pmap, 'degradationSet': degr_set,
                'degradationInputs': degr_inputs}


class RbaProcessingMap:

    def __init__(self, pmid, const_processing=None, comp_processing=None):
        self.id = pmid
        self.constant_processing = const_processing if type(const_processing) is dict else {}
        self.component_processings = comp_processing if type(comp_processing) is dict else {}

    @staticmethod
    def import_xml(processing_maps):

        data = {}
        for processing_map in processing_maps.findall('processingMap'):
            pmid = processing_map.attrib['id']
            rba_pmap = RbaProcessingMap(pmid)

            const_proc = processing_map.find('constantProcessing')
            if const_proc is not None:
                rba_pmap.constant_processing['reactants'] = get_species_refs_from_xml(
                    const_proc.find('listOfReactants'))
                rba_pmap.constant_processing['products'] = get_species_refs_from_xml(const_proc.find('listOfProducts'))

            comp_procs = processing_map.find('listOfComponentProcessings')
            if comp_procs is not None:
                for comp_proc in comp_procs.findall('componentProcessing'):
                    component = comp_proc.attrib['component']
                    cost = float(comp_proc.get('machineryCost', '0'))
                    reactants = get_species_refs_from_xml(comp_proc.find('listOfReactants'))
                    products = get_species_refs_from_xml(comp_proc.find('listOfProducts'))
                    rba_pmap.component_processings[component] = {'cost': cost, 'reactants': reactants,
                                                                 'products': products}
            data[pmid] = rba_pmap
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for pmid, group in df.groupby('processingMap'):
            rba_pmap = RbaProcessingMap(pmid)
            for _, row in group.iterrows():
                component = row['component']
                reactants = get_species_refs_from_str(row['reactants'])
                products = get_species_refs_from_str(row['products'])
                if component == 'constantProcessing':
                    rba_pmap.constant_processing['reactants'] = reactants
                    rba_pmap.constant_processing['products'] = products
                else:
                    cost = row['machineryCost']
                    rba_pmap.component_processings[component] = {'cost': cost, 'reactants': reactants,
                                                                 'products': products}
            data[pmid] = rba_pmap
        return data

    def export_xml(self):
        pmap = Element('processingMap', {'id': self.id})

        if len(self.constant_processing) > 0:
            const_proc = SubElement(pmap, 'constantProcessing')
            if len(self.constant_processing['reactants']) > 0:
                reactants = SubElement(const_proc, 'listOfReactants')
                for sid, stoic in self.constant_processing['reactants'].items():
                    SubElement(reactants, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
            if len(self.constant_processing['products']) > 0:
                products = SubElement(const_proc, 'listOfProducts')
                for sid, stoic in self.constant_processing['products'].items():
                    SubElement(products, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})

        if len(self.component_processings) > 0:
            comp_procs = SubElement(pmap, 'listOfComponentProcessings')
            for component, params in self.component_processings.items():
                attribs = {'component': component, 'machineryCost': str(params['cost'])}
                comp_proc = SubElement(comp_procs, 'componentProcessing', attribs)
                if len(params['reactants']) > 0:
                    reactants = SubElement(comp_proc, 'listOfReactants')
                    for sid, stoic in params['reactants'].items():
                        SubElement(reactants, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
                if len(params['products']) > 0:
                    products = SubElement(comp_proc, 'listOfProducts')
                    for sid, stoic in params['products'].items():
                        SubElement(products, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
        return pmap

    def to_df(self):
        pmid = self.id
        data = []
        if ((len(self.constant_processing.get('reactants', {})) > 0) or
                (len(self.constant_processing.get('reactants', {}))) > 0):
            cost = np.nan
            component = 'constantProcessing'
            reactants = '; '.join([f'species={species}, stoic={stoic}'
                                   for species, stoic in self.constant_processing['reactants'].items()])
            products = '; '.join([f'species={species}, stoic={stoic}'
                                  for species, stoic in self.constant_processing['products'].items()])
            data.append([pmid, component, cost, reactants, products])
        if len(self.component_processings) > 0:
            for component, comp_proc in self.component_processings.items():
                cost = comp_proc['cost']
                reactants = '; '.join([f'species={species}, stoic={stoic}'
                                       for species, stoic in comp_proc['reactants'].items()])
                products = '; '.join([f'species={species}, stoic={stoic}'
                                      for species, stoic in comp_proc['products'].items()])
                data.append([pmid, component, cost, reactants, products])
        df = pd.DataFrame(data, columns=['processingMap', 'component', 'cost', 'reactants', 'products'])
        return df

    def to_dict(self):
        const_proc = {}
        if 'reactants' in self.constant_processing:
            const_proc['reactants'] = '; '.join([f'species={species}, stoic={stoic}'
                                                 for species, stoic in self.constant_processing['reactants'].items()])
            const_proc['products'] = '; '.join([f'species={species}, stoic={stoic}'
                                                for species, stoic in self.constant_processing['products'].items()])

        comp_procs = {}
        if len(self.component_processings) > 0:
            for component, comp_proc in self.component_processings.items():
                mach_cost = comp_proc['cost']
                reactants = '; '.join([f'species={species}, stoic={stoic}'
                                       for species, stoic in comp_proc['reactants'].items()])
                products = '; '.join([f'species={species}, stoic={stoic}'
                                      for species, stoic in comp_proc['products'].items()])
                comp_procs[component] = {'machineryCost': mach_cost, 'reactants': reactants, 'products': products}

        return {'processing_map': self.id, 'constantProcessing': const_proc,
                'componentProcessings': comp_procs}

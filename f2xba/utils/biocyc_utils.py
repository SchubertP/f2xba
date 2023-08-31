"""Implementation of utility for bicocyc data parsing functions.

Peter Schubert, HHU Duesseldorf, May 2022
"""
import re
import numpy as np


def get_child_text(el_parent, child_tag):
    """retrieve child text data for a specified tag.

    There can be several elements with same tag,
    Text entries of all elements are retrieved and separated by '|'
    """
    text = ''
    if el_parent.find(child_tag) is not None:
        text = '| '.join([el_obj.text.strip() for el_obj in el_parent.findall(child_tag)])
    return strip_html_tags(text)


def to_float(string):
    try:
        return float(string)
    except ValueError:
        return np.nan


def to_int(string):
    try:
        return int(string)
    except ValueError:
        return 0


def strip_html_tags(string):
    return re.sub(r'<[^>]*>', '', string)


# get subtags of a specific component
def get_subtags(el_parent, child_tag):
    tags = []
    for el_child in el_parent.findall(child_tag):
        sub_tags = []
        for el_obj in el_child.findall('*'):
            sub_tags.append(el_obj.tag)
        tags.append('|'.join(sub_tags))
    return tags


def get_sub_obj_ids(el_parent, child_tag, obj_tag):
    obj_ids = []
    for el_child in el_parent.findall(child_tag):
        if el_child.find(obj_tag) is None:
            text = el_child.text
            if text is not None:
                obj_ids.append('(' + text.strip() + ')')
        else:
            for el_obj in el_child.findall(obj_tag):
                obj_ids.append(el_obj.get('frameid'))
    return obj_ids


def get_components(el_parent, component_type):
    """Retrieve components with stoichiometry

    :param el_parent: xml parent element from where to look for components
    :type el_parent: xml.etree.ElementTree.Element
    :param component_type: type of component to get information from
    :type  component_type: str
    :return: component ids with stoichiometry
    :rtype: dict (key: element id, value: stoichiometry as float)
    """
    components = {}
    for el_child in el_parent.findall('component'):
        ec_id = ''
        stoic = '1'
        for el_obj in el_child:
            if el_obj.tag == component_type:
                ec_id = el_obj.get('frameid')
            if el_obj.tag == 'coefficient':
                stoic = el_obj.text.strip()
        if ec_id != '':
            components[ec_id] = float(stoic)
    return components


# get resources
# type of resource (first 3 digits, compartment, coefficient)
def get_resources(el_parent, child_tag):
    resources = []
    for el_child in el_parent.findall(child_tag):
        ec_id = ''
        stoic = '1'
        compartment = ''
        for el_obj in el_child:
            if el_obj.tag in ['Protein', 'Compound', 'RNA']:
                ec_id = el_obj.get('frameid')
            if el_obj.tag == 'coefficient':
                stoic = el_obj.text.strip()
            if el_obj.tag == 'compartment':
                el_cco = el_obj.find('cco')
                compartment = el_cco.get('frameid')
        if ec_id != '':
            if compartment == '':
                resources.append(f'{ec_id}:{stoic}')
            else:
                resources.append(f'{ec_id}:{stoic}:{compartment}')
    return resources


def get_gene_products(el_parent, child_tag, product_type):
    """retrieve from XML tree a specified product type

    Note, there can be several protein or RNA products to a gene
    """
    products = []
    el_child = el_parent.find(child_tag)
    if el_child is not None:
        for el_obj in el_child:
            if el_obj.tag == product_type:
                product = el_obj.get('frameid')
                products.append(product)
    return products


def get_items(items_str, delim=';'):
    if type(items_str) != str:
        return []
    else:
        items = items_str.split(delim)
        for item in items:
            yield item.strip()


def extract_params(record):
    params = {}
    kv_pairs = record.split(',')
    for kv_pair in kv_pairs:
        if '=' in kv_pair:
            k, v = kv_pair.split('=')
            params[k.strip()] = v.strip()
    return params

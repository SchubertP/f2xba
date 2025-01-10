"""Implementation MiriamAnnotation class.


Peter Schubert, HHU Duesseldorf, December 2024
"""

from sbmlxdf.misc import record_generator
from collections import defaultdict


class MiriamAnnotation:

    def __init__(self, annot_str):
        self.references = {}
        if type(annot_str) is str:
            for annotation in record_generator(annot_str):
                fields = [item.strip() for item in annotation.split(',')]
                qualifier = fields[0]
                elements = defaultdict(set)
                for resource_uri in fields[1:]:
                    if '/' in resource_uri:
                        resource, data = resource_uri.rsplit('/', 1)
                    else:
                        resource, data = ('', resource_uri)
                    elements[resource].add(data)
                self.references[qualifier] = dict(elements)

    def get_annot_str(self):
        records = []
        for qualifier, qual_elements in self.references.items():
            qual_refs_elements = [qualifier]
            for resource, datas in qual_elements.items():
                for data in sorted(datas):
                    qual_refs_elements.append(f'{resource}/{data}')
            records.append(', '.join(qual_refs_elements))
        return '; '.join(records)

    def get_qualified_refs(self, qualifier, resource):
        if qualifier in self.references and resource in self.references[qualifier]:
            return sorted(self.references[qualifier][resource])
        else:
            return []

    def replace_refs(self, qualifier, resource, refs):
        if type(refs) is str:
            refs = [refs]
        if qualifier in self.references:
            self.references[qualifier][resource] = set(refs)
        else:
            self.references[qualifier] = {resource: set(refs)}

    def add_ref(self, qualifier, resource, ref):
        if qualifier in self.references:
            self.references[qualifier][resource].add(ref)
        else:
            self.references[qualifier] = {resource: {ref}}

    def delete(self, qualifier, resource=None, ref=None):
        if qualifier in self.references:
            if resource is None:
                del self.references[qualifier]
            else:
                if resource in self.references[qualifier]:
                    if ref is None:
                        del self.references[qualifier][resource]
                    else:
                        self.references[qualifier][resource].discard(ref)

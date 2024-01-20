"""Implementation of RbaTargetValue class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""


class RbaTargetValue:

    def __init__(self, value=None, lb=None, ub=None):
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub

    @staticmethod
    def get_target_value(value_type, p_name):
        """Return a RbaTargetValue object instantiated

        value_type can contain either 'value', 'upperBound' or 'lowerBound'
        allowed is also 'upperBound, lowerBound'

        :param value_type: value type
        :type value_type: str
        :param p_name: parameter name (function / aggretate)
        :type p_name: str
        :return: Target Value object
        :rType: RbaTargetValue
        """
        value = None
        lb = None
        ub = None
        if 'value' in value_type:
            value = p_name
        else:
            if 'upperBound' in value_type:
                ub = p_name
            if 'lowerBound' in value_type:
                lb = p_name
        return RbaTargetValue(value=value, lb=lb, ub=ub)

    @staticmethod
    def from_dict(value_dict):
        target_value = RbaTargetValue()
        if 'value' in value_dict:
            target_value.value = value_dict['value']
        else:
            target_value.lower_bound = value_dict.get('lowerBound')
            target_value.upper_bound = value_dict.get('upperBound')
        return target_value

    def to_dict(self):
        target_values = {}
        if self.lower_bound is not None:
            target_values['lowerBound'] = self.lower_bound
        if self.upper_bound is not None:
            target_values['upperBound'] = self.upper_bound
        if self.value is not None:
            target_values['value'] = self.value
        return target_values

    def to_str(self):
        return ', '.join([f'{key}={val}' for key, val in self.to_dict().items()])

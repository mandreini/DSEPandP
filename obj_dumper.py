__author__ = 'Matt'

from abc import ABCMeta
import pickle

class NullObject(object):
    pass

class MyObj:

    __metaclass__ = ABCMeta

    def __init__(self):
        self.a = 3
        self.b = 4
        self.c = 5

def obj_to_dict(obj):
    """
    Converts an object to a dictionary with key-values of the name of the methods/members as keys and
    the method/member themselves as the values. Requires classes to be an 'ABCmeta' __metaclass__ type.
    See: http://www.jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/
    for a brief turotial on ABCmeta classes.
    :param obj: Object to convert
    :return:
    """

    is_abcmeta = False

    if hasattr(obj, '__metaclass__'):
        if not obj.__metaclass__ == ABCMeta:
            not_abcmeta = True
    else:
        not_abcmeta = True

    if is_abcmeta:
        raise Exception("Object must be an ABCmeta class to comply with python new-style classes.")

    obj_dict = obj.__dict__
    obj_dict['__class__'] = obj.__class__

    return obj_dict

def dict_to_obj(obj_dict):
    """
    Converts a dictionary back into an object created from obj_to_dict method
    :param obj_dict: Dictionary corresponding to an object
    :return:
    """

    if not obj_dict.has_key('__class__'): raise Exception("This is not a dictionary created from the obj_to_dict procedure.")

    obj_class = obj_dict['__class__']
    new_obj = NullObject()
    new_obj.__class__ = obj_class

    for attr in obj_dict.keys():
        new_obj.__setattr__(attr, obj_dict[attr])

    return new_obj

def make_copy(obj):
    """
    Makes a copy (without using copy.copy) by calling obj_to_dict and then dict_to_obj
    :param obj: obj to copy
    :return: copy of obj
    """
    copy_dict = obj_to_dict(obj)
    copy_obj = dict_to_obj(copy_dict)
    return copy_obj



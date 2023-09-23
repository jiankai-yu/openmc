from collections.abc import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, c_size_t
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError, OpenMCError
from . import _dll, Nuclide
from .core import _FortranObjectWithID
from .error import _error_handler


__all__ = ['Surface', 'surfaces']

# Surface functions
_dll.openmc_extend_surfaces.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_surfaces.restype = c_int
_dll.openmc_extend_surfaces.errcheck = _error_handler
#
_dll.openmc_get_surface_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_surface_index.restype = c_int
_dll.openmc_get_surface_index.errcheck = _error_handler
#
_dll.openmc_surface_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_surface_get_id.restype = c_int
_dll.openmc_surface_get_id.errcheck = _error_handler
#
_dll.openmc_surface_get_coeff.argtypes = [c_int32, c_int32, POINTER(c_double)]
_dll.openmc_surface_get_coeff.restype = c_int
_dll.openmc_surface_get_coeff.errcheck = _error_handler
#
_dll.openmc_surface_get_name.argtypes = [c_int32, POINTER(c_char_p)]
_dll.openmc_surface_get_name.restype = c_int
_dll.openmc_surface_get_name.errcheck = _error_handler
#
_dll.openmc_surface_get_type.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_surface_get_type.restype = c_int
_dll.openmc_surface_get_type.errcheck = _error_handler
#
_dll.openmc_surface_set_id.argtypes = [c_int32, c_double, c_char_p]
_dll.openmc_surface_set_id.restype = c_int
_dll.openmc_surface_set_id.errcheck = _error_handler
#
_dll.openmc_surface_set_coeff.argtypes = [
    c_int32, c_int, c_double]
_dll.openmc_surface_set_coeff.restype = c_int
_dll.openmc_surface_set_coeff.errcheck = _error_handler
#
_dll.openmc_surface_set_name.argtypes = [c_int32, c_int32]
_dll.openmc_surface_set_name.restype = c_int
_dll.openmc_surface_set_name.errcheck = _error_handler
#
_dll.openmc_surface_set_type.argtypes = [c_int32, c_double]
_dll.openmc_surface_set_type.restype = c_int
_dll.openmc_surface_set_type.errcheck = _error_handler
#
_dll.n_surfaces.argtypes = []
_dll.n_surfaces.restype = c_size_t


class Surface(_FortranObjectWithID):
    """Surface stored internally.

    This class exposes a Surface that is stored internally in the OpenMC
    library. To obtain a view of a Surface with a given ID, use the
    :data:`openmc.lib.Surfaces` mapping.

    Parameters
    ----------
    uid : int or None
        Unique ID of the Surface
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view to an existing object is returned.
    index : int or None
         Index in the `Surfaces` array.

    Attributes
    ----------
    id : int
        ID of the Surface
    nuclides : list of str
        List of nuclides in the Surface
    densities : numpy.ndarray
        Array of densities in atom/b-cm
    name : str
        Name of the Surface
    temperature : float
        Temperature of the Surface in [K]
    volume : float
        Volume of the Surface in [cm^3]

    """
    __instances = WeakValueDictionary()

    def __new__(cls, uid=None, new=True, index=None):
        mapping = surfaces
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    uid = max(mapping, default=0) + 1
                else:
                    if uid in mapping:
                        raise AllocationError('A Surface with ID={} has already '
                                              'been allocated.'.format(uid))

                index = c_int32()
                _dll.openmc_extend_surfaces(1, index, None)
                index = index.value
            else:
                index = mapping[uid]._index
        elif index == -1:
            # Special value indicates void Surface
            return None

        if index not in cls.__instances:
            instance = super(Surface, cls).__new__(cls)
            instance._index = index
            if uid is not None:
                instance.id = uid
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        mat_id = c_int32()
        _dll.openmc_surface_get_id(self._index, mat_id)
        return mat_id.value

    @id.setter
    def id(self, mat_id):
        _dll.openmc_surface_set_id(self._index, mat_id)

    @property
    def name(self):
        name = c_char_p()
        _dll.openmc_surface_get_name(self._index, name)
        return name.value.decode()

    @name.setter
    def name(self, name):
        name_ptr = c_char_p(name.encode())
        _dll.openmc_surface_set_name(self._index, name_ptr)
    
    @property
    def type(self):
        type = c_char_p()
        _dll.openmc_surface_get_type(self._index, type)
        return type.value.decode()

    @name.setter
    def type(self, type):
        type_ptr = c_char_p(type.encode())
        _dll.openmc_surface_set_type(self._index, type_ptr)

    def get_coeff(self, n=1):
        """Get n-th coefficient of a Surface.

        Parameters
        ----------
        n : integral
            n-th coefficient of a surface 

        Returns
        -------
        float
            the requested n-th coefficient

        """
        coeff = c_double()
        _dll.openmc_surface_get_coeff(self._index, n, coeff)
        return coeff.value

    def set_coeff(self, n, coeff):
        """Set n-th coeffieicent of a Surface.

        Parameters
        ----------
        n : integral 
            n-th coefficient
        coeff : float
            coefficient of a surface

        """
        _dll.openmc_surface_set_coeff(self._index, n, coeff)


class _SurfaceMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_surface_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return Surface(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield Surface(index=i).id

    def __len__(self):
        return _dll.n_Surfaces()

    def __repr__(self):
        return repr(dict(self))

surfaces = _SurfaceMapping()

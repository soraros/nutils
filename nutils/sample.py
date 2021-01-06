# Copyright (c) 2014 Evalf
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

'''
The sample module defines the :class:`Sample` class, which represents a
collection of discrete points on a topology and is typically formed via
:func:`nutils.topology.Topology.sample`. Any function evaluation starts from
this sampling step, which drops element information and other topological
properties such as boundaries and groups, but retains point positions and
(optionally) integration weights. Evaluation is performed by subsequent calls
to :func:`Sample.integrate`, :func:`Sample.integral` or :func:`Sample.eval`.

Besides the location of points, :class:`Sample` also keeps track of point
connectivity through its :attr:`Sample.tri` and :attr:`Sample.hull`
properties, representing a (n-dimensional) triangulation of the interior and
boundary, respectively. Availability of these properties depends on the
selected sample points, and is typically used in combination with the "bezier"
set.
'''

from . import types, points, util, function, evaluable, parallel, numeric, matrix, transformseq, sparse
from .pointsseq import PointsSequence
import numpy, numbers, collections.abc, os, treelog as log, abc

graphviz = os.environ.get('NUTILS_GRAPHVIZ')

def argdict(arguments):
  if len(arguments) == 1 and 'arguments' in arguments and isinstance(arguments['arguments'], collections.abc.Mapping):
    arguments = arguments['arguments']
  return types.frozendict[types.strictstr,types.frozenarray](arguments)

class Sample(types.Singleton):
  '''Collection of points on a topology.

  The :class:`Sample` class represents a collection of discrete points on a
  topology and is typically formed via :func:`nutils.topology.Topology.sample`.
  Any function evaluation starts from this sampling step, which drops element
  information and other topological properties such as boundaries and groups,
  but retains point positions and (optionally) integration weights. Evaluation
  is performed by subsequent calls to :func:`integrate`, :func:`integral` or
  :func:`eval`.

  Besides the location of points, ``Sample`` also keeps track of point
  connectivity through its :attr:`tri` and :attr:`hull` properties,
  representing a (n-dimensional) triangulation of the interior and boundary,
  respectively. Availability of these properties depends on the selected sample
  points, and is typically used in combination with the "bezier" set.
  '''

  __slots__ = 'nelems', 'transforms', 'points', 'ndims'
  __cache__ = 'allcoords'

  @staticmethod
  @types.apply_annotations
  def new(transforms:types.tuple[transformseq.stricttransforms], points:types.strict[PointsSequence], index:types.tuple[types.frozenarray[int]]=None):
    '''Create a new :class:`Sample`.

    Parameters
    ----------
    transforms : :class:`tuple` or transformation chains
        List of transformation chains leading to local coordinate systems that
        contain points.
    points : :class:`~nutils.pointsseq.PointsSequence`
        Points sequence.
    index : :class:`tuple` of integer arrays, optional
        List of indices matching ``transforms``, defining the order on which
        points show up in the evaluation. If absent the indices will be strict
        increasing.
    '''

    if index is None:
      return _DefaultIndex(transforms, points)
    else:
      return _CustomIndex(transforms, points, index)

  def __init__(self, transforms, points):
    '''
    parameters
    ----------
    transforms : :class:`tuple` or transformation chains
        List of transformation chains leading to local coordinate systems that
        contain points.
    points : :class:`~nutils.pointsseq.PointsSequence`
        Points sequence.
    '''

    assert len(transforms) >= 1
    assert all(len(t) == len(points) for t in transforms)
    self.nelems = len(transforms[0])
    self.transforms = transforms
    self.points = points
    self.ndims = transforms[0].fromdims

  def __repr__(self):
    return '{}.{}<{}D, {} elems, {} points>'.format(type(self).__module__, type(self).__qualname__, self.ndims, self.nelems, self.npoints)

  @property
  def npoints(self):
    return self.points.npoints

  @property
  def index(self):
    return tuple(map(self.getindex, range(self.nelems)))

  @abc.abstractmethod
  def getindex(self, ielem):
    '''Return the indices of `Sample.points[ielem]` in results of `Sample.eval`.'''

    raise NotImplementedError

  def get_evaluable_indices(self, ielem):
    '''Return the evaluable indices for the given evaluable element index.

    Parameters
    ----------
    ielem : :class:`nutils.evaluable.Array`, ndim: 0, dtype: :class:`int`
        The element index.

    Returns
    -------
    indices : :class:`nutils.evaluable.Array`
        The indices of the points belonging to the given element as a 1D array.

    See Also
    --------
    :meth:`getindex` : the non-evaluable equivalent
    '''

    return evaluable.ElemwiseFromCallable(self.getindex, ielem, self.points.get_evaluable_coords(ielem).shape[:1], int)

  @util.positional_only
  @util.single_or_multiple
  @types.apply_annotations
  def integrate(self, funcs, arguments:argdict=...):
    '''Integrate functions.

    Args
    ----
    funcs : :class:`nutils.function.Array` object or :class:`tuple` thereof.
        The integrand(s).
    arguments : :class:`dict` (default: None)
        Optional arguments for function evaluation.
    '''

    datas = self.integrate_sparse(funcs, arguments)
    with log.iter.fraction('assembling', datas) as items:
      return [_convert(data, inplace=True) for data in items]

  @util.single_or_multiple
  @types.apply_annotations
  def integrate_sparse(self, funcs:types.tuple[function.asarray], arguments:types.frozendict[str,types.frozenarray]=None):
    '''Integrate functions into sparse data.

    Args
    ----
    funcs : :class:`nutils.function.Array` object or :class:`tuple` thereof.
        The integrand(s).
    arguments : :class:`dict` (default: None)
        Optional arguments for function evaluation.
    '''

    return eval_integrals_sparse(*map(self.integral, funcs), **(arguments or {}))

  def integral(self, func):
    '''Create Integral object for postponed integration.

    Args
    ----
    func : :class:`nutils.function.Array`
        Integrand.
    '''

    return _Integral(func, self)

  @util.positional_only
  @util.single_or_multiple
  @types.apply_annotations
  def eval(self, funcs, arguments:argdict=...):
    '''Evaluate function.

    Args
    ----
    funcs : :class:`nutils.function.Array` object or :class:`tuple` thereof.
        The integrand(s).
    arguments : :class:`dict` (default: None)
        Optional arguments for function evaluation.
    '''

    datas = self.eval_sparse(funcs, arguments)
    with log.iter.fraction('assembling', datas) as items:
      return tuple(map(sparse.toarray, items))

  @util.positional_only
  @util.single_or_multiple
  @types.apply_annotations
  def eval_sparse(self, funcs:types.tuple[function.asarray], arguments:types.frozendict[str,types.frozenarray]=None):
    '''Evaluate function.

    Args
    ----
    funcs : :class:`nutils.function.Array` object or :class:`tuple` thereof.
        The integrand(s).
    arguments : :class:`dict` (default: None)
        Optional arguments for function evaluation.
    '''

    return eval_integrals_sparse(*map(self.__rmatmul__, funcs), **(arguments or {}))

  def __rmatmul__(self, func: function.Array) -> function.Array:
    if not isinstance(func, function.Array):
      return NotImplemented
    return _AtSample(func, self)

  @property
  def allcoords(self):
    coords = numpy.empty([self.npoints, self.ndims])
    for ielem, points in enumerate(self.points):
      coords[self.getindex(ielem)] = points.coords
    return types.frozenarray(coords, copy=False)

  def basis(self):
    '''Basis-like function that for every point in the sample evaluates to the
    unit vector corresponding to its index.'''

    index = function.transforms_index(self.transforms[0])
    coords = function.transforms_coords(self.transforms[0], self.ndims)
    I = function.Elemwise(self.index, index, dtype=int)
    B = function.Sampled(coords, expect=function.take(self.allcoords, I, axis=0))
    return function.inflate(B, indices=I, length=self.npoints, axis=0)

  def asfunction(self, array):
    '''Convert sampled data to evaluable array.

    Using the result of :func:`Sample.eval`, create a
    :class:`nutils.function.Sampled` array that upon evaluation recovers the
    original function in the set of points matching the original sampling.

    >>> from nutils import mesh
    >>> domain, geom = mesh.rectilinear([1,2])
    >>> gauss = domain.sample('gauss', 2)
    >>> data = gauss.eval(geom)
    >>> sampled = gauss.asfunction(data)
    >>> domain.integrate(sampled, degree=2)
    array([ 1.,  2.])

    Args
    ----
    array :
        The sampled data.
    '''

    return function.matmat(self.basis(), array)

  @property
  def tri(self):
    '''Triangulation of interior.

    A two-dimensional integer array with ``ndims+1`` columns, of which every
    row defines a simplex by mapping vertices into the list of points.
    '''

    return numpy.concatenate([self.getindex(ielem).take(points.tri) for ielem, points in enumerate(self.points)])

  @property
  def hull(self):
    '''Triangulation of the exterior hull.

    A two-dimensional integer array with ``ndims`` columns, of which every row
    defines a simplex by mapping vertices into the list of points. Note that
    the hull often does contain internal element boundaries as the
    triangulations originating from separate elements are disconnected.
    '''

    return numpy.concatenate([self.getindex(ielem).take(points.hull) for ielem, points in enumerate(self.points)])

  def subset(self, mask):
    '''Reduce the number of points.

    Simple selection mechanism that returns a reduced Sample based on a
    selection mask. Points that are marked True will still be part of the new
    subset; points marked False may be dropped but this is not guaranteed. The
    point order of the original Sample is preserved.

    Args
    ----
    mask : :class:`bool` array.
        Boolean mask that selects all points that should remain. The resulting
        Sample may contain more points than this, but not less.

    Returns
    -------
    subset : :class:`Sample`
    '''

    selection = types.frozenarray([ielem for ielem in range(self.nelems) if mask[self.getindex(ielem)].any()])
    transforms = tuple(transform[selection] for transform in self.transforms)
    return Sample.new(transforms, self.points.take(selection))

strictsample = types.strict[Sample]

class _DefaultIndex(Sample):

  __slots__ = ()
  __cache__ = 'offsets'

  @property
  def offsets(self):
    return numpy.cumsum([0]+[p.npoints for p in self.points])

  def getindex(self, ielem):
    return numpy.arange(self.offsets[ielem], self.offsets[ielem+1])

  @property
  def tri(self):
    return self.points.tri

  @property
  def hull(self):
    return self.points.hull

  def get_evaluable_indices(self, ielem):
    npoints = self.points.get_evaluable_coords(ielem).shape[0]
    offsets = evaluable._SizesToOffsets(evaluable.loop_concatenate(evaluable.InsertAxis(npoints, 1), ielem, self.nelems))
    return evaluable.Range(npoints, evaluable.get(offsets, 0, ielem))

class _CustomIndex(Sample):

  __slots__ = '_index'

  def __init__(self, transforms, points, index):
    self._index = index
    if len(index) != len(points):
      raise ValueError('expected an `index` with {} items but got {}'.format(len(points), len(index)))
    if not all(len(i) == p.npoints for i, p in zip(self._index, points)):
      raise ValueError('lengths of indices does not match number of points per element')
    super().__init__(transforms, points)

  @property
  def index(self):
    return self._index

  def getindex(self, ielem):
    return self._index[ielem]

@types.apply_annotations
def eval_integrals(*integrals: types.tuple, **arguments:argdict):
  '''Evaluate integrals.

  Evaluate one or several postponed integrals. By evaluating them
  simultaneously, rather than using :meth:`nutils.function.Array.eval` on each
  integral individually, integrations will be grouped per Sample and jointly
  executed, potentially increasing efficiency.

  Args
  ----
  integrals : :class:`tuple` of integrals
      Integrals to be evaluated.
  arguments : :class:`dict` (default: None)
      Optional arguments for function evaluation.

  Returns
  -------
  results : :class:`tuple` of arrays and/or :class:`nutils.matrix.Matrix` objects.
  '''

  with log.iter.fraction('assembling', eval_integrals_sparse(*integrals, **arguments)) as retvals:
    return [_convert(retval, inplace=True) for retval in retvals]

@types.apply_annotations
def eval_integrals_sparse(*integrals: types.tuple, **arguments: argdict):
  '''Evaluate integrals into sparse data.

  Evaluate one or several postponed integrals. By evaluating them
  simultaneously, rather than using :meth:`nutils.function.Array.eval` on each
  integral individually, integrations will be grouped per Sample and jointly
  executed, potentially increasing efficiency.

  Args
  ----
  integrals : :class:`tuple` of integrals
      Integrals to be evaluated.
  arguments : :class:`dict` (default: None)
      Optional arguments for function evaluation.

  Returns
  -------
  results : :class:`tuple` of arrays and/or :class:`nutils.matrix.Matrix` objects.
  '''

  # it = integrals[0]
  # subit = tuple(it.func.func.funcs)[1]
  # # print(tuple(subit)[1])
  # print(subit.dtype)
  # print('here>>>>>')
  # print(subit.asciitree())
  # print(it.asciitree())
  # print(it.simplified.asciitree())
  # it.simplified.eval(**arguments)
  integrals = tuple(integral.as_evaluable_array().assparse for integral in integrals)
  with evaluable.Tuple(tuple(integrals)).optimized_for_numpy.session(graphviz=graphviz) as eval:
    return eval(**arguments)

def _convert(data, inplace=False):
  '''Convert a two-dimensional sparse object to an appropriate object.

  The return type is determined based on dimension: a zero-dimensional object
  becomes a scalar, a one-dimensional object a (dense) Numpy vector, a
  two-dimensional object a Nutils matrix, and any higher dimensional object a
  deduplicated and pruned sparse object.
  '''

  ndim = sparse.ndim(data)
  return sparse.toarray(data) if ndim < 2 \
    else matrix.fromsparse(data, inplace=inplace) if ndim == 2 \
    else sparse.prune(sparse.dedup(data, inplace=inplace), inplace=True)

class _Integral(function.Array):

  def __init__(self, integrand: function.Array, sample: Sample) -> None:
    self._integrand = integrand
    self._sample = sample
    # print('good', integrand.dtype)
    super().__init__(shape=integrand.shape, dtype=float if integrand.dtype in (bool, int) else integrand.dtype)

  def lower(self, *, transform_chains=None, coordinates=None, **kwargs) -> evaluable.Array:
    if transform_chains or coordinates:
      raise ValueError('nested integrals or samples are not yet supported')
    ielem = evaluable.Argument('_ielem', (), dtype=int)
    transform_chains = tuple(evaluable.TransformChainFromSequence(t, ielem) for t in self._sample.transforms)
    coordinates = (self._sample.points.get_evaluable_coords(ielem),) * len(self._sample.transforms)
    integrand = self._integrand.lower(transform_chains=transform_chains, coordinates=coordinates, **kwargs)
    return evaluable.LoopSum(evaluable.dot(evaluable.appendaxes(self._sample.points.get_evaluable_weights(ielem), integrand.shape[1:]), integrand, 0), ielem, self._sample.nelems)

class _AtSample(function.Array):

  def __init__(self, func: function.Array, sample: Sample) -> None:
    self._func = func
    self._sample = sample
    super().__init__(shape=(sample.points.npoints, *func.shape), dtype=func.dtype)

  def lower(self, *, transform_chains=None, coordinates=None, **kwargs) -> evaluable.Array:
    if transform_chains or coordinates:
      raise ValueError('nested integrals or samples are not yet supported')
    ielem = evaluable.Argument('_ielem', (), dtype=int)
    transform_chains = tuple(evaluable.TransformChainFromSequence(t, ielem) for t in self._sample.transforms)
    coordinates = (self._sample.points.get_evaluable_coords(ielem),) * len(self._sample.transforms)
    func = self._func.lower(transform_chains=transform_chains, coordinates=coordinates)
    indices = self._sample.get_evaluable_indices(ielem)
    inflated = evaluable.Transpose.from_end(evaluable.Inflate(evaluable.Transpose.to_end(func, 0), indices, self._sample.npoints), 0)
    return evaluable.LoopSum(inflated, ielem, self._sample.nelems)

# vim:sw=2:sts=2:et

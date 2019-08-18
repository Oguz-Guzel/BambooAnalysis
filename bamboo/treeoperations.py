"""
Object representation of operations on TTree branches

The aim is to provide provide sufficiently flexible and complete foundation
for the development of efficient histogram-filling programs
through the use of python wrappers (see e.g. treeproxies).
"""

from itertools import chain, repeat, combinations, count, tee
from contextlib import contextmanager
import logging
logger = logging.getLogger(__name__)

from collections import defaultdict
_perfCnt = defaultdict(int)
def _printPerfCnt():
    logger.info("Counters of hits per method in TupleOp")
    for nm, cnt in sorted(_perfCnt.items(), key=(lambda elm : elm[1]), reverse=True):
        logger.info(" - {0:>10d} {1}".format(cnt, nm))

def simpletrace(func):
    def wrapped_fun(self, *args, **kwargs):
        _perfCnt[".".join((self.__class__.__name__, func.__name__))] += 1
        return func(self, *args, **kwargs)
    wrapped_fun.__name__ = func.__name__
    wrapped_fun.__doc__ = func.__doc__
    return wrapped_fun

class TupleOpCache:
    __slots__ = ("hash", "repr")
    def __init__(self):
        self.hash = None
        self.repr = None
    def __bool__(self):
        return self.hash is not None or self.repr is not None

def fromopcache(func):
    key = func.__name__.strip("_")
    def wrapped_prop(self):
        perfNm = ".".join((self.__class__.__name__, key))
        if getattr(self._cache, key) is None:
            _perfCnt["{0}_w".format(perfNm)] += 1
            setattr(self._cache, key, func(self))
        else:
            _perfCnt["{0}_c".format(perfNm)] += 1
        return getattr(self._cache, key)
    wrapped_prop.__name__ = func.__name__
    wrapped_prop.__doc__ = func.__doc__
    return wrapped_prop

class TupleOp:
    """ Interface & base class for operations on leafs and resulting objects / values

    Instances should be defined once, and assumed immutable by all observers
    (they should only ever be modified just after construction, preferably by the owner).
    Since a value-based hash (and repr) is cached on first use, violating this rule
    might lead to serious bugs. In case of doubt the clone method can be used to
    obtain an independent copy.
    Subclasses should define a result property and clone, _repr, _eq, and optionally deps methods.
    """
    __slots__ = ("_cache") ## this means all deriving classes need to define __slots__ (for performance)
    def __init__(self):
        self._cache = TupleOpCache()
    def clone(self, memo=None):
        """ Create an independent copy (with empty repr/hash cache) of the (sub)expression """
        if memo is None: ## top-level, construct the dictionary
            memo = dict()
        if id(self) in memo:
            return memo[id(self)]
        else:
            cp = self._clone(memo)
            memo[id(self)] = cp
            return cp
    def _clone(self, memo): ## simple version, call clone of attributes without worrying about memo
        """ Implementation of clone - to be overridden by all subclasses (memo is dealt with by clone, so simply construct, calling .clone(memo=memo) on TupleOp attributes """
        return self.__class__()
    def deps(self, defCache=None, select=(lambda x : True), includeLocal=False):
        """ Dependent TupleOps iterator """
        yield from []
    @property
    def result(self):
        """ Proxy to the result of this (sub)expression """
        pass
    ## subclasses should define at least _clone, _repr, and _eq (value-based)
    @fromopcache
    def __repr__(self):
        """ String representation (used for hash, and lazily cached) """
        return self._repr()
    def _repr(self):
        """ __repr__ implementation - to be overridden by all subclasses (caching is in top-level __repr__) """
        return "TupleOp()"
    @fromopcache
    def __hash__(self):
        """ Value-based hash (lazily cached) """
        return hash(self.__repr__())
    @simpletrace
    def __eq__(self, other):
        """ Identity or value-based equality comparison (same object and unequal should be fast) """
        # _eq may end up being quite expensive, but should almost never be called
        return id(self) == id(other) or ( self.__hash__() == hash(other) and self.__class__ == other.__class__ and self._eq(other) )
    def _eq(self, other):
        """ value-based __eq__ implementation - to be overridden by all subclasses (protects against hash collisions; hash and class are checked to be equal already) """
        return True

## implementations are split out, see treeproxies
class TupleBaseProxy:
    """
    Interface & base class for proxies
    """
    def __init__(self, typeName, parent=None):
        self._typeName = typeName
        self._parent = parent
    @property
    def op(self):
        if self._parent is None:
            raise ValueError("Cannot get operation for {0!r}, abstract base class / empty parent".format(self))
        return self._parent

class CppStrRedir:
    """ Expression cache interface. Default implementation: no caching """
    def __init__(self):
        self._iFun = 0
    def __call__(self, arg):
        return arg.get_cppStr(defCache=self)
    def symbol(self, decl):
        """
        Define (or get) a new C++ symbol for the declaration

        decl should contain the code, with <<name>> where the name should go.  Returns the unique name
        """
        print("WARNING: should add defined symbol for '{0}' but that's not supported".format(decl))
    def _getColName(self, op):
        return None

cppNoRedir = CppStrRedir()

class ForwardingOp(TupleOp):
    """ Transparent wrapper (base for marking parts of the tree, e.g. things with systematic variations) """
    __slots__ = ("wrapped",)
    def __init__(self, wrapped):
        super(ForwardingOp, self).__init__()
        self.wrapped = wrapped
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.wrapped.clone(memo))
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from self.wrapped.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        return self.wrapped.result
    def _repr(self):
        return "{0}({1!r})".format(self.__class__.__name__, self.wrapped)
    def _eq(self, other):
        return self.wrapped == other.wrapped
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return self.wrapped.get_cppStr(defCache=defCache)

SizeType = "std::size_t"

class Const(TupleOp):
    """ Hard-coded number (or expression) """
    __slots__ = ("typeName", "value")
    def __init__(self, typeName, value):
        super(Const, self).__init__()
        self.typeName = typeName
        self.value = value
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, self.value)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.value)
    def _eq(self, other):
        return self.typeName == other.typeName and self.value == other.value
    # backends
    @simpletrace
    def get_cppStr(self, defCache=None):
        try:
            if abs(self.value) == float("inf"):
                return "std::numeric_limits<{0}>::{mnmx}()".format(self.typeName, mnmx=("min" if self.value < 0. else "max"))
        except:
            pass
        return str(self.value) ## should maybe be type-aware...

class GetColumn(TupleOp):
    """ Get a column value """
    __slots__ = ("typeName", "name")
    def __init__(self, typeName, name):
        super(GetColumn, self).__init__()
        self.typeName = typeName
        self.name = name
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, self.name)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.name)
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name
    @simpletrace
    def get_cppStr(self, defCache=None):
        return self.name

class GetArrayColumn(TupleOp):
    """ Get the number from a leaf """
    __slots__ = ("typeName", "name", "length")
    def __init__(self, typeName, name, length):
        super(GetArrayColumn, self).__init__()
        self.typeName = typeName
        self.name = name
        self.length = length
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, self.name, self.length.clone(memo=memo))
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if select(self.length):
            yield self.length
        yield from self.length.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self, makeProxy(SizeType, self.length))
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.typeName, self.name, self.length)
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name and self.length == other.length
    @simpletrace
    def get_cppStr(self, defCache=None):
        return self.name

## helper
def adaptArg(arg, typeHint=None):
    if isinstance(arg, TupleBaseProxy):
        return arg.op
    elif isinstance(arg, TupleOp):
        return arg
    elif typeHint is not None:
        if str(arg) == arg: ## string, needs quote
            return Const(typeHint, '"{}"'.format(arg))
        else:
            return Const(typeHint, arg)
    else:
        raise ValueError("Should get some kind of type hint")

mathOpFuns_cppStr = {
      "add"      : lambda cppStr,*args : "( {0} )".format(" + ".join(cppStr(arg) for arg in args))
    , "multiply" : lambda cppStr,*args : "( {0} )".format(" * ".join(cppStr(arg) for arg in args))
    , "subtract" : lambda cppStr,a1,a2 : "( {0} - {1} )".format(cppStr(a1), cppStr(a2))
    , "divide"   : lambda cppStr,a1,a2 : "( {0} / {1} )".format(cppStr(a1), cppStr(a2))
    , "floatdiv"  : lambda cppStr,a1,a2 : "( 1.*{0} / {1} )".format(cppStr(a1), cppStr(a2))
    #
    , "lt" : lambda cppStr,a1,a2 : "( {0} <  {1} )".format(cppStr(a1), cppStr(a2))
    , "le" : lambda cppStr,a1,a2 : "( {0} <= {1} )".format(cppStr(a1), cppStr(a2))
    , "eq" : lambda cppStr,a1,a2 : "( {0} == {1} )".format(cppStr(a1), cppStr(a2))
    , "ne" : lambda cppStr,a1,a2 : "( {0} != {1} )".format(cppStr(a1), cppStr(a2))
    , "gt" : lambda cppStr,a1,a2 : "( {0} >  {1} )".format(cppStr(a1), cppStr(a2))
    , "ge" : lambda cppStr,a1,a2 : "( {0} >= {1} )".format(cppStr(a1), cppStr(a2))
    , "and" : lambda cppStr,*args : "( {0} )".format(" && ".join(cppStr(a) for a in args))
    , "or"  : lambda cppStr,*args : "( {0} )".format(" || ".join(cppStr(a) for a in args))
    , "not" : lambda cppStr,a : "( ! {0} )".format(cppStr(a))
    , "band" : lambda cppStr,*args : "( {0} )".format(" & ".join(cppStr(a) for a in args))
    , "bor"  : lambda cppStr,*args : "( {0} )".format(" | ".join(cppStr(a) for a in args))
    , "bxor"  : lambda cppStr,*args : "( {0} )".format(" ^ ".join(cppStr(a) for a in args))
    , "bnot" : lambda cppStr,a : "( ~ {0} )".format(cppStr(a))
    #
    , "abs"   : lambda cppStr,arg : "std::abs( {0} )".format(cppStr(arg))
    , "sqrt"  : lambda cppStr,arg : "std::sqrt( {0} )".format(cppStr(arg))
    , "pow"   : lambda cppStr,a1,a2 : "std::pow( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "exp"   : lambda cppStr,arg : "std::exp( {0} )".format(cppStr(arg))
    , "log"   : lambda cppStr,arg : "std::log( {0} )".format(cppStr(arg))
    , "log10" : lambda cppStr,arg : "std::log10( {0} )".format(cppStr(arg))
    , "max"   : lambda cppStr,a1,a2 : "std::max( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "min"   : lambda cppStr,a1,a2 : "std::min( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    #
    , "switch" : lambda cppStr,test,trueBr,falseBr : "( {0} ) ? ( {1} ) : ( {2} )".format(cppStr(test), cppStr(trueBr), cppStr(falseBr))
    }

class MathOp(TupleOp):
    """ Mathematical function N->1, e.g. sin, abs, ( lambda x, y : x*y ) """
    __slots__ = ("outType", "op", "args")
    def __init__(self, op, *args, **kwargs):
        super(MathOp, self).__init__()
        self.outType = kwargs.pop("outType", "Double_t")
        assert len(kwargs) == 0
        self.op = op
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.op, *(a.clone(memo=memo) for a in self.args), outType=self.outType)
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.outType, self)
    def _repr(self):
        return "{0}({1}, {2}, outType={3!r})".format(self.__class__.__name__, self.op, ", ".join(repr(arg) for arg in self.args), self.outType)
    def _eq(self, other):
        return self.outType == other.outType and self.op == other.op and self.args == other.args
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return mathOpFuns_cppStr[self.op](defCache, *self.args)

class GetItem(TupleOp):
    """ Get item from array (from function call or from array leaf) """
    __slots__ = ("arg", "typeName", "_index")
    def __init__(self, arg, valueType, index, indexType=SizeType):
        super(GetItem, self).__init__()
        self.arg = adaptArg(arg)
        self.typeName = valueType
        self._index = adaptArg(index, typeHint=indexType)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.arg.clone(memo=memo), self.typeName, self._index.clone(memo=memo))
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for arg in (self.arg, self._index):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def index(self):
        from .treeproxies import makeProxy
        return makeProxy(SizeType, self._index)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.arg, self.typeName, self._index)
    def _eq(self, other):
        return self.arg == other.arg and self.typeName == other.typeName and self._index == other._index
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}[{1}]".format(defCache(self.arg), defCache(self._index))

class Construct(TupleOp):
    __slots__ = ("typeName", "args")
    def __init__(self, typeName, args):
        super(Construct, self).__init__()
        self.typeName = typeName
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, tuple(a.clone(memo=memo) for a in self.args))
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2})".format(self.__class__.__name__, self.typeName, ", ".join(repr(a) for a in self.args))
    def _eq(self, other):
        return self.typeName == other.typeName and self.args == other.args
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}{{{1}}}".format(self.typeName, ", ".join(defCache(a) for a in self.args))

def guessReturnType(mp):
    if hasattr(mp, "func_doc") and hasattr(mp, "func_name"):
        toks = list(mp.func_doc.split())
        ## left and right strip const * and &
        while toks[-1].rstrip("&") in ("", "const", "static"):
            toks = toks[:-1]
        while toks[0].rstrip("&") in ("", "const", "static"):
            toks = toks[1:]
        while any(tok.endswith("unsigned") for tok in toks):
            iU = next(i for i,tok in enumerate(toks) if tok.endswith("unsigned"))
            toks[iU] = " ".join((toks[iU], toks[iU+1]))
            del toks[iU+1]
        if len(toks) == 2:
            return toks[0].rstrip("&")
        else:
            nOpen = 0
            i = 0
            while i < len(toks) and ( i == 0 or nOpen != 0 ):
                nOpen += ( toks[i].count("<") - toks[i].count(">") )
                i += 1
            return " ".join(toks[:i]).rstrip("&")
    else:
        return "Float_t"

class CallMethod(TupleOp):
    """
    Call a method
    """
    __slots__ = ("name", "args", "_retType")
    def __init__(self, name, args, returnType=None, getFromRoot=True):
        super(CallMethod, self).__init__()
        self.name = name ## NOTE can only be a hardcoded string this way
        self.args = tuple(adaptArg(arg) for arg in args)
        self._retType = returnType if returnType else CallMethod._initReturnType(name, getFromRoot=getFromRoot)
    @staticmethod
    def _initReturnType(name, getFromRoot=True):
        mp = None
        if getFromRoot:
            try:
                from cppyy import gbl
                if "::" in name:
                    res = gbl
                    for tok in name.split("::"):
                        res = getattr(res, tok)
                    if res != gbl:
                        mp = res
                else:
                    mp = getattr(gbl, name)
            except Exception as ex:
                logger.error("Exception in getting method pointer {0}: {1}".format(name, ex))
        return guessReturnType(mp)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.name, tuple(a.clone(memo=memo) for a in self.args), returnType=self._retType)
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self._retType, self)
    def _repr(self):
        return "{0}({1!r}, ({2}), returnType={3!r})".format(self.__class__.__name__, self.name, ", ".join(repr(arg) for arg in self.args), self._retType)
    def _eq(self, other):
        return self.name == other.name and self._retType == other._retType and self.args == other.args
    # backends
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        if not defCache.shouldDefine(self):
            return "{0}({1})".format(self.name, ", ".join(defCache(arg) for arg in self.args))
        else: ## go through a symbol
            depList = _collectDeps(self.args, [], defCache=defCache)
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = "{name}({args})\n".format(name=self.name, args=", ".join(defCache(arg) for arg in self.args))
            funName = defCache.symbol(expr, resultType=self.result._typeName, args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class CallMemberMethod(TupleOp):
    """ Call a member method """
    __slots__ = ("this", "name", "args", "_retType")
    def __init__(self, this, name, args, returnType=None):
        super(CallMemberMethod, self).__init__()
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        self.args = tuple(adaptArg(arg) for arg in args)
        self._retType = returnType if returnType else guessReturnType(getattr(this._typ, self.name))
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.this.clone(memo=memo), self.name, tuple(a.clone(memo=memo) for a in self.args), returnType=self._retType)
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in chain((self.this,), self.args):
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self._retType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, ({3}), returnType={4!r})".format(self.__class__.__name__, self.this, self.name, ", ".join(repr(arg) for arg in self.args), self._retType)
    def _eq(self, other):
        return self.this == other.this and self.name == other.name and self._retType == other._retType and self.args == other.args
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}({2})".format(defCache(self.this), self.name, ", ".join(defCache(arg) for arg in self.args))

class GetDataMember(TupleOp):
    """ Get a data member """
    __slots__ = ("this", "name")
    def __init__(self, this, name):
        super(GetDataMember, self).__init__()
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.this.clone(memo=memo), self.name)
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            if select(self.this):
                yield self.this
            yield from self.this.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        if not self.name.startswith("_"):
            try:
                protoTp = self.this.result._typ
                proto = protoTp() ## should *in principle* work for most ROOT objects
                att = getattr(proto, self.name)
                tpNm = type(att).__name__
                if protoTp.__name__.startswith("pair<") and self.name in ("first", "second"):
                    tpNms = tuple(tok.strip() for tok in protoTp.__name__[5:-1].split(","))
                    return makeProxy((tpNms[0] if self.name == "first" else tpNms[1]), self)
                return makeProxy(tpNm, self)
            except Exception as e:
                print("Problem getting type of data member {0} of {1!r}".format(self.name, self.this), e)
        return makeProxy("void", self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.this, self.name)
    def _eq(self, other):
        return self.this == other.this and self.name == other.name
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}".format(defCache(self.this), self.name)

class ExtVar(TupleOp):
    """ Externally-defined variable (used by name) """
    __slots__ = ("typeName", "name")
    def __init__(self, typeName, name):
        super(ExtVar, self).__init__()
        self.typeName = typeName
        self.name = name
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, self.name)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.name)
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name
    @simpletrace
    def get_cppStr(self, defCache=None):
        return self.name

class DefinedVar(TupleOp):
    """ Defined variable (used by name), first use will trigger definition """
    __slots__ = ("typeName", "definition", "_nameHint")
    def __init__(self, typeName, definition, nameHint=None):
        super(DefinedVar, self).__init__()
        self.typeName = typeName
        self.definition = definition
        self._nameHint = nameHint
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, self.definition, nameHint=self._nameHint)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, nameHint={3!r})".format(self.__class__.__name__, self.typeName, self.definition, self._nameHint)
    def _eq(self, other):
        return self.typeName == other.typeName and self.definition == other.definition and self._nameHint == other._nameHint
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return defCache.symbol(self.definition, nameHint=self._nameHint)

class InitList(TupleOp):
    """ Initializer list """
    __slots__ = ("typeName", "elms")
    def __init__(self, typeName, elms, elmType=None):
        super(InitList, self).__init__()
        self.typeName = typeName
        self.elms = tuple(adaptArg(e, typeHint=elmType) for e in elms)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeName, tuple(elm.clone(memo=memo) for elm in self.elms))
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for elm in self.elms:
                if select(elm):
                    yield elm
                yield from elm.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}<{1}>({2})".format(self.__class__.__name__, self.typeName, ", ".join(repr(elm) for elm in self.elms))
    def _eq(self, other):
        return self.typeName == other.typeName and self.elms == other.elms
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        return "{{ {0} }}".format(", ".join(defCache(elm) for elm in self.elms))

class LocalVariablePlaceholder(TupleOp):
    """ Placeholder type for a local variable connected to an index (first step in a specific-to-general strategy) """
    __slots__ = ("typeHint", "_parent", "i")
    def __init__(self, typeHint, parent=None, i=None):
        super(LocalVariablePlaceholder, self).__init__()
        self.typeHint = typeHint
        self._parent = parent
        self.i = i ## FIXME this one is set **late** - watch out with what we call
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.typeHint, parent=self._parent, i=self.i)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeHint, self)
    @property
    def name(self):
        if self.i is None:
            raise RuntimeError("Using LocalVariablePlaceholder before giving it an index")
        return "i{0:d}".format(self.i)
    @simpletrace
    def get_cppStr(self, defCache=None):
        return self.name
    def _repr(self):
        return "{0}({1!r}, i={2!r})".format(self.__class__.__name__, self.typeHint, self.i)
    def _eq(self, other):
        ## NOTE this breaks the infinite recursion, but may not be 100% safe
        ## what should save the nested case is that the repr(parent) will be different for different levels of nesting
        ## since all LVP's are supposed to have an index, confusion between cases where they are combined differently should be eliminated as well
        return self.typeHint == other.typeHint and repr(self._parent) == repr(other._parent) and self.i == other.i

def collectNodes(expr, select=(lambda nd : True)):
    # simple helper
    if select(expr):
        yield expr
    yield from expr.deps(select=select, includeLocal=True)

def _collectDeps(exprs, ownLocal, defCache=cppNoRedir):
    ## first pass (will trigger definitions, if necessary)
    exprs1, exprs2 = tee(exprs, 2)
    for dep in chain.from_iterable(expr.deps(defCache=defCache, select=lambda op : defCache.shouldDefine(op)) for expr in exprs1):
        cn = defCache(dep)
        if not cn:
            logger.warning("Probably a problem in triggering definition for {0}".format(dep))
    return set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn)
                or defCache.shouldDefine(op) or ( isinstance(op, LocalVariablePlaceholder) and op not in ownLocal )
                ))
            for expr in exprs2))

def _convertFunArgs(deps, defCache=cppNoRedir):
    capDeclCall = []
    for ld in deps:
        if isinstance(ld, GetArrayColumn):
            capDeclCall.append((
                "&{0}".format(ld.name),
                "const ROOT::VecOps::RVec<{0}>& {1}".format(ld.typeName, ld.name),
                ld.name))
        elif isinstance(ld, GetColumn):
            capDeclCall.append((
                "&{0}".format(ld.name),
                "const {0}& {1}".format(ld.typeName, ld.name),
                ld.name))
        elif isinstance(ld, LocalVariablePlaceholder):
            if not ld.name:
                print("ERROR: no name for local {0}".format(ld))
            capDeclCall.append((
                ld.name,
                "{0} {1}".format(ld.typeHint, ld.name),
                ld.name))
        elif defCache.shouldDefine(ld):
            nm = defCache._getColName(ld)
            if not nm:
                print("ERROR: no column name for {0}".format(ld))
            if not any("&{0}".format(nm) == icap for icap,idecl,icall in capDeclCall):
                capDeclCall.append((
                    "&{0}".format(nm),
                    "const {0}& {1}".format(ld.result._typeName, nm),
                    nm))
            else:
                print("WARNING: dependency {0} is there twice".format(nm))
        else:
            raise AssertionError("Dependency with unknown type: {0}".format(ld))
    return zip(*sorted(capDeclCall, key=(lambda elm : elm[1]))) ## sort by declaration (alphabetic for the same type)

def _normFunArgs(expr, args, argNames):
    newExpr = expr
    newArgs = args
    for i,argN in sorted(list(enumerate(argNames)), key=(lambda elm : len(elm[1])), reverse=True):
        newName = "myArg{0:d}".format(i, argN)
        assert sum(1 for ia in newArgs if argN in ia) == 1 ## should be in one and only one argument
        newArgs = [ (ia.replace(argN, newName) if argN in ia else ia) for ia in newArgs ]
        newExpr = newExpr.replace(argN, newName)
    return newExpr, newArgs

class Select(TupleOp):
    """ Define a selection on a range """
    __slots__ = ("rng", "predExpr", "_i")
    def __init__(self, rng, predExpr, idx):
        super(Select, self).__init__()
        self.rng = rng ## proxy
        self.predExpr = predExpr
        self._i = idx
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.rng.clone(memo=memo), self.predExpr.clone(memo=memo), self._i.clone(memo=memo))
    @staticmethod
    def fromRngFun(rng, pred):
        """ Factory method from a range and predicate (callable) """
        idx = LocalVariablePlaceholder(SizeType)
        predExpr = adaptArg(pred(rng._base[idx.result]))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(predExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Select(rng, predExpr, idx)
        idx._parent = res
        return res
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (adaptArg(self.rng), self.predExpr):
                if select(arg):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(SizeType), itemType=SizeType)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.predExpr, self._i)
    def _eq(self, other):
        return self.rng == other.rng and self.predExpr == other.predExpr and self._i == other._i
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.predExpr), (self._i,), defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::select({idxs},\n    [{captures}] ( {i} ) {{ return {predExpr}; }})".format(
                idxs=defCache(self.rng._idxs.op),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                predExpr=defCache(self.predExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, args_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType="ROOT::VecOps::RVec<{0}>".format(SizeType), args=", ".join(args_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Sort(TupleOp):
    """ Sort a range (ascendingly) by the value of a function on each element """
    __slots__ = ("rng", "funExpr", "_i")
    def __init__(self, rng, funExpr, idx):
        super(Sort, self).__init__()
        self.rng = rng ## PROXY
        self.funExpr = funExpr
        self._i = idx
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.rng.clone(memo=memo), self.funExpr.clone(memo=memo), self._i.clone(memo=memo))
    @staticmethod
    def fromRngFun(rng, fun):
        idx = LocalVariablePlaceholder(SizeType)
        funExpr = adaptArg(fun(rng._base[idx.result]))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(funExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Sort(rng, funExpr, idx)
        idx._parent = res
        return res
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (adaptArg(self.rng), self.funExpr):
                if select(arg):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(SizeType), itemType=SizeType)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.funExpr, self._i)
    def _eq(self, other):
        return self.rng == other.rng and self.funExpr == other.funExpr and self._i == other._i
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.funExpr), (self._i,), defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::sort({idxs},\n    [{captures}] ( {i} ) {{ return {funExpr}; }})".format(
                idxs=defCache(self.rng._idxs.op),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                funExpr=defCache(self.funExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            funName = defCache.symbol(expr, resultType="ROOT::VecOps::RVec<{0}>".format(SizeType), args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Map(TupleOp):
    """ Create a list of derived values for a collection (mostly useful for storing on skims) """
    __slots__ = ("rng", "funExpr", "_i", "typeName")
    def __init__(self, rng, funExpr, idx, typeName):
        super(Map, self).__init__()
        self.rng = rng ## PROXY
        self.funExpr = funExpr
        self._i = idx
        self.typeName = typeName
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.rng.clone(memo=memo), self.funExpr.clone(memo=memo), self._i.clone(memo=memo), self.typeName)
    @staticmethod
    def fromRngFun(rng, fun, typeName=None):
        idx = LocalVariablePlaceholder(SizeType)
        val = fun(rng._base[idx.result])
        funExpr = adaptArg(val)
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(funExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Map(rng, funExpr, idx, typeName=(typeName if typeName is not None else val._typeName))
        idx._parent = res
        return res
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (adaptArg(self.rng), self.funExpr):
                if select(arg):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(self.typeName), itemType=self.typeName)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self.rng, self.funExpr, self._i, self.typeName)
    def _eq(self, other):
        return self.rng == other.rng and self.funExpr == other.funExpr and self._i == other._i and self.typeName == other.typeName
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.funExpr), (self._i,), defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::map<{valueType}>({idxs},\n    [{captures}] ( {i} ) {{ return {funExpr}; }})".format(
                valueType=self.typeName,
                idxs=defCache(self.rng._idxs.op),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                funExpr=defCache(self.funExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            funName = defCache.symbol(expr, resultType="ROOT::VecOps::RVec<{0}>".format(self.typeName), args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Next(TupleOp):
    """ Define a search (first matching item, for a version that processes the whole range see Reduce) """
    __slots__ = ("rng", "predExpr", "_i")
    def __init__(self, rng, predExpr, idx):
        super(Next, self).__init__()
        self.rng = rng ## PROXY
        self.predExpr = predExpr
        self._i = idx
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.rng.clone(memo=memo), self.predExpr.clone(memo=memo), self._i.clone(memo=memo))
    @staticmethod
    def fromRngFun(rng, pred): ## FIXME you are here
        idx = LocalVariablePlaceholder(SizeType)
        predExpr = adaptArg(pred(rng._base[idx.result]))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(predExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Next(rng, predExpr, idx)
        idx._parent = res
        return res
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (adaptArg(self.rng), self.predExpr):
                if select(arg):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        return self.rng._base[self]
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.predExpr, self._i)
    def _eq(self, other):
        return self.rng == other.rng and self.predExpr == other.predExpr and self._i == other._i
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.predExpr), (self._i,), defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::next({idxs},\n     [{captures}] ( {i} ) {{ return {predexpr}; }}, -1)".format(
                idxs=defCache(self.rng._idxs.op),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                predexpr=defCache(self.predExpr),
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            funName = defCache.symbol(expr, resultType=SizeType, args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Reduce(TupleOp):
    """ Reduce a range to a value (could be a transformation, index...) """
    def __init__(self, rng, resultType, start, accuExpr, idx, prevRes):
        self.rng = rng ## PROXY
        self.resultType = resultType
        self.start = start
        self.accuExpr = accuExpr
        self._i = idx
        self._prevRes = prevRes
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.rng.clone(memo=memo), self.resultType, self.start.clone(memo=memo), self.accuExpr.clone(memo=memo), self._i.clone(memo=memo), self._prevRes.clone(memo=memo))
    @staticmethod
    def fromRngFun(rng, start, accuFun):
        resultType = start._typeName
        idx = LocalVariablePlaceholder(SizeType)
        prevRes = LocalVariablePlaceholder(resultType, i=-1)
        accuExpr = adaptArg(accuFun(prevRes.result, rng._base[idx.result]))
        maxLVIdx = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(accuExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))
        idx.i = maxLVIdx+1
        prevRes.i = maxLVIdx+2

        res = Reduce(rng, resultType, adaptArg(start), accuExpr, idx, prevRes)
        idx._parent = res
        prevRes._parent = res
        return res
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.start, self.accuExpr):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in (self._i, self._prevRes):
                        yield dp
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.resultType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r}, {5!r}, {6!r})".format(self.__class__.__name__, self.rng, self.resultType, self.start, self.accuExpr, self._i, self_prevRes)
    def _eq(self, other):
        return self.rng == other.rng and self.resultType == other.resultType and self.start == other.start and self.accuExpr == other.accuExpr and self._i == other._i and self._prevRes == other._prevRes
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.start, self.accuExpr), (self._i, self._prevRes), defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::reduce({idxs}, {start},\n     [{captures}] ( {prevRes}, {i} ) {{ return {accuexpr}; }})".format(
                idxs=defCache(self.rng._idxs.op),
                start=defCache(self.start),
                captures=", ".join(captures),
                prevRes="{0} {1}".format(self._prevRes.typeHint, self._prevRes.name),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                accuexpr=defCache(self.accuExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            funName = defCache.symbol(expr, resultType=self.resultType, args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Combine(TupleOp):
    __slots__ = ("ranges", "candPredExpr", "_i")
    def __init__(self, ranges, candPredExpr, idx):
        super(Combine, self).__init__()
        self.ranges = ranges ## (PROXY,)
        self.candPredExpr = candPredExpr
        self._i = idx
    @property
    def n(self):
        return len(self.ranges)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(tuple(rng.clone(memo=memo) for rng in self.ranges), self.candPredExpr.clone(memo=memo), tuple(i.clone(memo=memo) for i in self._i))
    @staticmethod
    def fromRngFun(num, ranges, candPredFun, sameIdxPred=lambda i1,i2: i1 < i2):
        ranges = ranges if len(ranges) > 1 else tuple(repeat(ranges[0], num))
        idx = tuple(LocalVariablePlaceholder(SizeType, i=-1-i) for i in range(num))
        from . import treefunctions as op
        areDiff = op.AND(*(sameIdxPred(ia.result, ib.result)
                for ((ia, ra), (ib, rb)) in combinations(zip(idx, ranges), 2)
                if ra._base == rb._base))
        candPred = candPredFun(*( rng._base[iidx.result] for rng,iidx in zip(ranges, idx)))
        if len(areDiff.op.args) > 0:
            candPredExpr = adaptArg(op.AND(areDiff, candPred))
        else:
            candPredExpr = adaptArg(candPred)
        maxLVIdx = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(candPredExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))
        for i,ilvp in enumerate(idx):
            ilvp.i = maxLVIdx+1+i
        res = Combine(ranges, candPredExpr, idx)
        for ilvp in idx:
            ilvp._parent = res
        return res
    @property
    def resultType(self):
        return "ROOT::VecOps::RVec<rdfhelpers::Combination<{0:d}>>".format(self.n)
    @simpletrace
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in chain(self.ranges, [self.candPredExpr]):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import CombinationListProxy, makeProxy
        return CombinationListProxy(self, makeProxy(self.resultType, self))
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.ranges, self.candPredExpr, self._i)
    def _eq(self, other):
        return self.ranges == other.ranges and self.candPredExpr == other.candPredExpr and self._i == other._i
    @simpletrace
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(chain(self.ranges, [self.candPredExpr]), self._i, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = ("rdfhelpers::combine{num:d}(\n"
            "     [{captures}] ( {predIdxArgs} ) {{ return {predExpr}; }},\n"
            "     {ranges})").format(
                num=self.n,
                captures=", ".join(captures),
                predIdxArgs=", ".join("{0} {1}".format(i.typeHint, i.name) for i in self._i),
                predExpr = defCache(self.candPredExpr),
                ranges=", ".join(defCache(rng._idxs.op) for rng in self.ranges)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            funName = defCache.symbol(expr, resultType=self.resultType, args=", ".join(paramDecl))
            return "{0}({1})".format(funName, ", ".join(paramCall))

## FIXME to be implemented
class PsuedoRandom(TupleOp):
    """ Pseudorandom number (integer or float) within range """
    def __init__(self, xMin, xMax, seed, isIntegral=False):
        super(PseudoRandom, self).__init__()
        self.xMin = xMin
        self.xMax = xMax
        self.seed = seed
        self.isIntegral = isIntegral
    @property
    def resultType(self):
        return "Int_" if self.isIntegral else "Float_t"
    ## deps from xMin, xMax and seed
    ## seed can be event-based or object-based, depending?
    ## TODO implement C++ side as well

class OpWithSyst(ForwardingOp):
    """ Interface and base class for nodes that can change the systematic variation of something they wrap """
    def __init__(self, wrapped, systName, variations=None):
        super(OpWithSyst, self).__init__(wrapped)
        self.systName = systName
        self.variations = variations
    def changeVariation(self, newVar):
        pass ## main interface method
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.wrapped, self.systName, self.variations)
    def _eq(self, other):
        return super(OpWithSyst, self)._eq(other) and self.systName == other.systName and self.variations == other.variations

class ScaleFactorWithSystOp(OpWithSyst):
    """ Scalefactor (ILeptonScaleFactor::get() call), to be modified with Up/Down variations (these are cached) """
    def __init__(self, wrapped, systName, variations=None):
        super(ScaleFactorWithSystOp, self).__init__(wrapped, systName, variations=(variations if variations else [ "{0}{1}".format(systName, vard) for vard in ["up", "down"] ]))
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.wrapped.clone(memo=memo), self.systName, variations=self.variations)
    def changeVariation(self, newVariation):
        """ Assumed to be called on a fresh copy - *will* change the underlying value """
        if self._cache: # validate this assumption
            raise RuntimeError("Cannot change variation of an expression that is already frozen")
        if newVariation not in self.variations:
            raise ValueError("Invalid variation: {0}".format(newVariation))
        newVariation = (newVariation[len(self.systName):] if newVariation.startswith(self.systName) else newVariation).capitalize() ## translate to name in C++
        if self.wrapped.args[-1].name == "Nominal" and newVariation != self.wrapped.args[-1].name:
            self.wrapped.args[-1].name = newVariation

class SystModifiedCollectionOp(OpWithSyst):
    """ modifiedcollections 'at' call, to be modified to get another collection """
    def __init__(self, wrapped, name, variations):
        super(SystModifiedCollectionOp, self).__init__(wrapped, name, variations=variations)
    @simpletrace
    def _clone(self, memo):
        return self.__class__(self.wrapped.clone(memo=memo), self.systName, list(self.variations))
    def changeVariation(self, newCollection):
        """ Assumed to be called on a fresh copy - *will* change the underlying value """
        if self._cache: # validate this assumption
            raise RuntimeError("Cannot change variation of an expression that is already frozen")
        if newCollection not in self.variations:
            raise ValueError("Invalid collection: {0}".format(newCollection))
        if self.wrapped.args[0].value == '"nominal"' and newCollection != self.wrapped.args[0].value.strip('"'):
            self.wrapped.args[0].value = '"{0}"'.format(newCollection)

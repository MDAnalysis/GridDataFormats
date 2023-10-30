Search.setIndex({"docnames": ["gridData/basic", "gridData/core", "gridData/formats", "gridData/formats/OpenDX", "gridData/formats/gOpenMol", "gridData/formats/mrc", "gridData/overview", "index", "installation"], "filenames": ["gridData/basic.rst", "gridData/core.rst", "gridData/formats.rst", "gridData/formats/OpenDX.rst", "gridData/formats/gOpenMol.rst", "gridData/formats/mrc.rst", "gridData/overview.rst", "index.rst", "installation.rst"], "titles": ["Basic use", "Core functionality for storing n-D grids \u2014 <code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">gridData.core</span></code>", "Formats", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">OpenDX</span></code> \u2014 routines to read and write simple OpenDX files", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">gOpenMol</span></code> \u2014 the gOpenMol plt format", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">mrc</span></code> \u2014 the MRC/CCP4 volumetric data format", "Handling grids of data \u2014 <code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">gridData</span></code>", "GridDataFormats: Handling volumetric data in Python", "Installation"], "terms": {"In": [0, 1, 3, 6], "most": [0, 3], "case": [0, 1, 4], "onli": [0, 1, 2, 3, 4, 5, 6], "one": [0, 1, 3, 4, 6], "class": [0, 2, 6, 7], "i": [0, 1, 2, 3, 4, 5, 6, 7, 8], "import": [0, 3, 6], "grid": [0, 2, 3, 5, 7], "so": [0, 1, 2, 3], "we": [0, 3, 6, 8], "just": [0, 1, 3], "thi": [0, 1, 3, 4, 5, 6, 8], "right": 0, "awai": 0, "from": [0, 1, 4, 5, 6, 7, 8], "griddata": [0, 3, 4, 5, 7], "opendx": [0, 1, 2, 4], "file": [0, 1, 5, 7], "g": [0, 1, 3, 4, 6], "dx": [0, 1, 2, 6], "see": [0, 1, 3, 4, 6, 7], "also": [0, 1, 3, 4, 6, 7, 8], "read": [0, 1, 2, 4, 5, 7, 8], "more": [0, 3, 4], "inform": [0, 3, 4], "especi": [0, 1, 6], "when": [0, 1, 3], "work": [0, 1, 3, 4, 5, 7], "visual": [0, 1, 2, 3], "program": [0, 3, 4, 6], "pymol": [0, 2, 3, 6], "vmd": [0, 2, 3, 6], "chimera": [0, 2, 3], "gopenmol": [0, 2], "plt": [0, 2], "output": [0, 1, 3, 6], "numpi": [0, 1, 4, 5, 6], "histogramdd": [0, 1, 3, 4, 5, 6], "r": 0, "random": 0, "randn": 0, "100": [0, 4], "3": [0, 1, 3, 4, 5, 7], "h": 0, "edg": [0, 1, 3, 4, 5, 6], "np": 0, "bin": [0, 1, 3], "5": [0, 1, 4], "8": [0, 4], "4": [0, 1, 3, 4], "For": [0, 1, 3, 6, 8], "other": [0, 1, 3, 4, 5, 6, 8], "wai": [0, 3, 6], "doc": [0, 3, 6], "some": [0, 1, 3, 6], "format": [0, 1, 3, 7], "support": [0, 1, 3, 5, 6, 7, 8], "detail": [0, 3, 6], "core": [0, 3], "export": [0, 1, 3, 6], "method": [0, 1, 3, 4], "The": [0, 1, 2, 3, 4, 5, 6, 7, 8], "can": [0, 1, 2, 3, 4, 5, 6, 7, 8], "specifi": [0, 1, 3], "explicitli": 0, "pkl": 0, "file_format": [0, 1, 3], "pickl": [0, 1, 2], "mai": [0, 1], "take": [0, 1, 3], "addit": [0, 3, 6], "specif": [0, 1, 3], "keyword": [0, 1, 3], "which": [0, 1, 3, 6], "ar": [0, 1, 2, 3, 4, 5, 6, 7, 8], "document": [0, 3, 4, 6, 8], "separ": [0, 3], "assum": [0, 1], "ha": [0, 3, 8], "were": [0, 2], "gener": [0, 1, 3, 7], "same": [0, 1, 6], "posit": [0, 1, 3, 6], "store": [0, 3, 6], "A": [0, 1, 2, 4, 6, 7], "b": [0, 1, 3], "first": [0, 1, 3, 4], "object": [0, 1, 6], "c": [0, 5], "result": [0, 1], "ani": [0, 1, 3, 4, 8], "capabl": 0, "viewer": [0, 2], "later": [0, 8], "again": [0, 6], "interpol": [0, 1, 6], "cubic": [0, 1], "spline": [0, 1], "twice": 0, "sampl": 0, "a2": 0, "resample_factor": [0, 1], "2": [0, 1, 3, 4, 5, 7], "downsampl": 0, "half": 0, "each": [0, 1, 3, 4, 6], "dimens": [0, 1, 3, 4, 6], "ahalf": 0, "0": [0, 1, 3, 4, 5, 6, 7], "anoth": [0, 1, 6], "a_on_b": 0, "even": 0, "simpler": [0, 1], "region": [0, 6], "valu": [0, 1, 3, 4], "did": 0, "occur": 0, "origin": [0, 1, 3, 4, 5, 6], "particular": [0, 1], "": [0, 1, 3], "lowest": 0, "wa": 0, "probabl": 0, "produc": [0, 1, 3, 4], "where": [0, 5], "chang": [0, 1, 3], "abruptli": 0, "modul": [1, 3, 4, 6, 7], "contain": [1, 3, 6, 7, 8], "independ": [1, 4], "data": [1, 2, 3], "act": [1, 6], "univers": [1, 3, 6], "constructor": [1, 3, 4, 6], "kwarg": [1, 3], "construct": [1, 3], "filenam": [1, 3, 4, 5, 6], "desir": [1, 6], "make": [1, 6, 7], "an": [1, 3, 6], "empti": 1, "load": [1, 3, 4, 5], "popul": [1, 4, 5], "none": [1, 3, 4, 5], "delta": [1, 3, 4, 5, 6], "metadata": 1, "interpolation_spline_ord": 1, "sourc": [1, 3, 4, 5, 7, 8], "multidimension": 1, "space": [1, 3, 5, 6], "us": [1, 2, 3, 4, 5, 7], "arithmet": [1, 6], "calcul": 1, "like": [1, 3, 4], "arrai": [1, 2, 4, 5, 6], "thei": 1, "compat": [1, 3], "e": [1, 3, 4, 5, 6], "have": [1, 3, 4, 6, 8], "shape": [1, 3, 4, 5], "length": [1, 3], "order": [1, 3, 5], "resampl": [1, 6], "common": [1, 3], "attribut": [1, 5], "hold": [1, 3], "standard": [1, 2, 4, 6], "directli": [1, 2, 6], "manipul": [1, 6], "number": [1, 2, 3, 4, 7], "molecular": [1, 2, 4], "volum": 1, "densiti": [1, 3, 4, 5, 6], "written": [1, 2, 3, 4, 6, 7], "out": [1, 4, 6], "differ": [1, 4], "paramet": [1, 3, 5], "ndarrai": [1, 3, 5, 6], "str": [1, 3, 5], "option": [1, 3, 5], "build": [1, 8], "either": [1, 3, 4], "histogram": [1, 3, 4, 5], "nd": [1, 3, 4, 5], "list": [1, 3, 7], "lower": [1, 3], "upper": 1, "along": [1, 3], "ax": [1, 3, 6], "cartesian": [1, 3, 6], "coordin": [1, 3, 4, 5, 6], "center": [1, 3], "index": [1, 3, 4, 7], "x": [1, 2, 4, 5], "cell": [1, 3, 4, 5, 6], "1": [1, 3, 4, 5, 7], "rectangular": [1, 3, 6], "dict": [1, 3, 4], "user": [1, 2], "defin": [1, 3, 4, 6], "dictionari": [1, 3], "arbitrari": [1, 3, 6], "kei": 1, "pair": 1, "associ": [1, 3], "doe": 1, "touch": 1, "save": [1, 6], "int": [1, 3, 5], "default": [1, 3], "name": [1, 3, 8], "necessari": [1, 8], "autodetect": 1, "fail": [1, 4, 8], "normal": 1, "guess": [1, 4], "extens": [1, 2, 3], "rais": [1, 3, 5], "typeerror": [1, 4], "If": [1, 3, 4, 5], "variou": 1, "input": [1, 3, 5], "do": [1, 3, 4, 8], "agre": 1, "valueerror": [1, 3, 5], "requir": [1, 3], "provid": [1, 3, 4], "argument": [1, 3], "suppli": [1, 4], "notimplementederror": 1, "triclin": 1, "non": 1, "orthorhomb": [1, 5], "box": 1, "note": [1, 3, 5, 6], "1d": 1, "ndim": 1, "repres": [1, 3, 4, 5, 6], "high": 1, "dimension": [1, 3, 5, 6], "real": [1, 3], "axi": [1, 3, 5], "convent": 1, "griddataformat": 1, "correspond": [1, 3, 5, 6], "compon": [1, 3], "y": [1, 4, 5], "z": [1, 4, 5], "type": [1, 3, 4, 5], "voxels": 1, "system": [1, 3, 4, 5, 6], "describ": [1, 3, 4, 7], "becaus": [1, 3, 4], "boundari": 1, "between": [1, 4], "all": [1, 3, 4, 6, 8], "last": [1, 4], "regular": [1, 3, 4, 5, 6, 7], "indic": 1, "midpoint": 1, "annot": 1, "content": 1, "It": [1, 3, 4, 6, 8], "togeth": [1, 6], "exampl": [1, 3, 4], "creat": [1, 3, 6], "principl": 1, "practic": 1, "mani": [1, 2, 3], "three": 1, "almost": 1, "test": [1, 4], "special": 1, "alwai": [1, 3, 4, 5], "3d": [1, 4, 5], "might": 1, "python": [1, 2, 3, 6, 8], "version": [1, 3, 5, 7, 8], "new": [1, 2, 3, 5, 6], "7": [1, 4, 5], "ccp4": [1, 2], "now": [1, 3], "mrc": [1, 2], "anymor": 1, "deprec": 1, "buggi": 1, "return": [1, 3, 4, 5], "iter": [1, 3], "ndindex": 1, "check_compat": 1, "check": [1, 3], "oper": [1, 6], "scalar": [1, 3], "float": [1, 3, 4], "default_format": 1, "typequot": [1, 3], "given": 1, "deduc": 1, "suffix": 1, "although": 1, "preced": 1, "implement": [1, 2, 3, 4], "restor": 1, "than": [1, 3], "set": [1, 3], "doubl": [1, 3], "By": [1, 3], "determin": [1, 3], "dtype": [1, 3], "typic": [1, 3], "charact": 1, "quot": [1, 3], "string": [1, 3, 6], "custom": 1, "parser": [1, 3], "namd": 1, "gridforc": 1, "backend": 1, "mdff": 1, "expect": 1, "appeas": 1, "them": [1, 6, 8], "properti": [1, 4, 5, 6], "over": [1, 6], "allow": [1, 3, 6], "obtain": 1, "x1": 1, "x2": 1, "y1": 1, "y2": 1, "z1": 1, "z2": 1, "f": 1, "comput": [1, 5], "onc": [1, 8], "cach": 1, "better": 1, "perform": 1, "whenev": 1, "modifi": [1, 4], "recomput": 1, "unknown": [1, 3], "interpolation_cv": 1, "todo": 1, "usag": 1, "xx": 1, "yy": 1, "zz": 1, "mgrid": 1, "40": 1, "75": 1, "96": 1, "150": 1, "20": 1, "50": [1, 4], "ff": 1, "possibl": [1, 4], "would": 1, "appear": 1, "neg": 1, "intern": 1, "scipi": [1, 8], "ndimag": 1, "map_coordin": 1, "mode": 1, "constant": [1, 6], "wherebi": 1, "outsid": 1, "fill": [1, 3], "beyond": 1, "minimum": 1, "6": [1, 3, 4, 7], "rather": 1, "nearest": 1, "elimin": 1, "extrud": 1, "choos": 1, "accept": 1, "spline_filt": 1, "call": [1, 3], "complet": [1, 3], "reset": 1, "base": 1, "current": [1, 3, 4, 5], "befor": 1, "tupl": 1, "instanc": [1, 3, 4, 5, 6], "As": 1, "conveni": [1, 3, 6], "othergrid": 1, "taken": [1, 5, 6], "factor": 1, "scale": 1, "n_i": 1, "must": [1, 3, 4], "cannot": 1, "fewer": 1, "size": [1, 3, 4, 5], "previou": 1, "alter": 1, "rang": 1, "being": 1, "creep": 1, "steadili": 1, "inward": 1, "recalcul": 1, "extent": 1, "everi": 1, "add": [1, 3, 6, 8], "regener": 1, "depend": [1, 3, 8], "therefor": 1, "guarante": [1, 3], "sai": 1, "altern": 1, "portabl": [1, 3], "ndmeshgrid": 1, "arr": 1, "mesh": [1, 6], "entri": 1, "fed": 1, "point": [1, 4, 6], "span": 1, "http": [1, 3, 4, 7], "stackoverflow": 1, "com": [1, 7], "question": [1, 7], "1827489": 1, "meshgrid": 1, "fix": [1, 3], "limit": [2, 7], "commonli": [2, 7], "particularli": 2, "suitabl": 2, "interfac": [2, 3], "tool": [2, 6], "ad": [2, 8], "difficult": 2, "contribut": [2, 7], "reader": [2, 3, 4, 6], "writer": 2, "easili": [2, 8], "integr": 2, "send": 2, "pull": [2, 7], "request": [2, 7], "packag": [2, 6, 7, 8], "extend": 2, "wide": 2, "understood": [2, 6], "suffici": [2, 3], "applic": 2, "encount": 2, "far": 2, "henc": 2, "moment": 2, "small": 2, "write": [2, 4, 6], "remark": 2, "subset": [2, 3], "routin": 2, "simpl": [2, 4, 6, 7], "dxinitobject": [2, 3], "dxparseerror": [2, 3], "dxparser": [2, 3], "dxparsernotoken": [2, 3], "dxclass": [2, 3], "field": [2, 3], "gridconnect": [2, 3], "gridposit": [2, 3], "volumetr": 2, "multi": 3, "free": 3, "softwar": [3, 5], "www": [3, 4], "org": [3, 4], "primit": 3, "n": [3, 6], "appendix": 3, "explor": 3, "nativ": 3, "you": [3, 4, 8], "want": [3, 4], "your": [3, 4, 7, 8], "top": 3, "level": 3, "below": 3, "realli": 3, "need": [3, 8], "should": [3, 4, 8], "tell": 3, "about": 3, "what": [3, 8], "howev": 3, "full": 3, "veri": [3, 4], "tri": [3, 4], "sometim": 3, "help": 3, "apb": 3, "seventh": 3, "signific": 3, "figur": [3, 4], "reflect": 3, "increas": 3, "precis": 3, "section": 3, "35": 3, "sinc": 3, "releas": [3, 7], "close": 3, "approxim": 3, "often": 3, "float64": 3, "forc": 3, "dx_type": 3, "avail": [3, 6, 7], "for_pymol": 3, "abl": [3, 4], "suggest": 3, "parlanc": 3, "unit": [3, 5], "follow": 3, "time": 3, "off": 3, "diagon": [3, 5], "element": 3, "ones": 3, "width": 3, "eg": 3, "angstrom": 3, "connect": 3, "d": [3, 4], "One": 3, "exist": 3, "initi": 3, "ident": 3, "account": 3, "manner": 3, "link": 3, "spec": 3, "sdsc": 3, "edu": 3, "html": 3, "page": [3, 4], "usrgu068": 3, "htm": 3, "hdredf": 3, "dead": 3, "am": 3, "archiv": [3, 4], "copi": [3, 7], "internet": [3, 4], "classtyp": 3, "classid": 3, "storag": 3, "variabl": [3, 4], "arg": 3, "turn": 3, "except": 3, "pars": 3, "error": 3, "brain": 3, "baroqu": 3, "scan": 3, "line": [3, 4, 6], "id": 3, "accord": 3, "setup": 3, "dxfield_object": 3, "p": 3, "bulk": 3, "rewritten": 3, "includ": 3, "found": [3, 4], "layout": 3, "ignor": 3, "remov": 3, "apply_pars": 3, "appli": 3, "token": 3, "stream": [3, 5], "dxfield": 3, "hierarchi": 3, "examin": 3, "distinguish": 3, "comment": 3, "basic": 3, "idea": 3, "state": 3, "machin": [3, 4], "There": 3, "activ": 3, "main": 3, "loop": [3, 4], "yet": 3, "set_pars": 3, "parsernam": 3, "use_pars": 3, "buffer": 3, "exhaust": 3, "ndformat": 3, "repetit": 3, "self": 3, "deriv": 3, "optstr": 3, "fals": 3, "pack": 3, "array_lik": 3, "cast": 3, "closest": 3, "under": [3, 7], "correct": [3, 8], "byte": 3, "uint8": 3, "float32": 3, "int32": 3, "short": 3, "int16": 3, "sign": 3, "int8": 3, "unsign": 3, "uint32": 3, "uint16": 3, "convers": 3, "round": 3, "trip": 3, "int64": 3, "miss": [3, 8], "np_type": 3, "float16": 3, "uint64": 3, "dxtype": 3, "valid": 3, "whole": [3, 8], "instanti": 3, "subclass": 3, "saniti": 3, "individu": 3, "becom": 3, "prefix": 3, "avoid": 3, "newlin": 3, "least": 3, "those": 3, "belong": 3, "dump": 3, "could": [3, 6], "keep": 3, "its": 3, "referenc": 3, "uniqu": 3, "gridpoint": 3, "dxobj": 3, "add_com": 3, "dxfile": 3, "discard": 3, "replac": 3, "sorted_compon": 3, "header": [3, 4, 5], "via": [3, 8], "statement": 3, "plugin": 3, "tuplet": [3, 4], "centr": [3, 4, 5], "dxd": [3, 4], "orthonorm": [3, 4, 5], "endia": 4, "own": 4, "csc": 4, "fi": 4, "english": 4, "g0penmol": 4, "develop": [4, 5, 7], "plt_format": 4, "phtml": 4, "access": 4, "through": [4, 7], "web": 4, "20061011125817": 4, "copyright": 4, "2005": 4, "septemb": 4, "23": 4, "2003": 4, "09": 4, "18": 4, "plot": 4, "orbit": 4, "electron": 4, "sever": 4, "unformat": 4, "pltfile": 4, "util": 4, "directori": 4, "extern": 4, "text": 4, "descript": 4, "elsewher": 4, "manual": 4, "hardwar": 4, "platform": 4, "littl": 4, "big": 4, "endian": 4, "pleas": [4, 7], "observ": 4, "fortran": 4, "pure": [4, 8], "record": [4, 5], "while": 4, "well": 4, "mean": 4, "integ": [4, 5], "rank": [4, 5], "surfac": 4, "vss": 4, "probe": 4, "200": 4, "gaussian": 4, "94": 4, "98": 4, "201": 4, "jaguar": 4, "202": 4, "gamess": 4, "203": 4, "autodock": 4, "204": 4, "delphi": 4, "insight": 4, "205": 4, "reserv": 4, "come": 4, "openmol": 4, "direct": [4, 5], "zmin": 4, "zmax": 4, "ymin": 4, "9": 4, "ymax": 4, "10": [4, 5], "xmin": 4, "11": 4, "xmax": 4, "12": 4, "run": [4, 8], "inner": 4, "nx": 4, "ny": 4, "nz": 4, "few": [4, 8], "look": 4, "65": 4, "300000e": 4, "001": 4, "200000e": 4, "625609e": 4, "644741e": 4, "663923e": 4, "683115e": 4, "702274e": 4, "721340e": 4, "740280e": 4, "759018e": 4, "777478e": 4, "795639e": 4, "813387e": 4, "830635e": 4, "zdim": 4, "ydim": 4, "xdim": 4, "min": 4, "max": 4, "per": 4, "held": 4, "mrc2014": 5, "mrcfile": 5, "librari": 5, "burnley2017": 5, "refer": 5, "burnlei": 5, "t": 5, "palmer": 5, "winn": 5, "m": 5, "2017": 5, "recent": 5, "ccp": 5, "em": 5, "suit": 5, "acta": 5, "cryst": 5, "d73": 5, "469": 5, "477": 5, "doi": 5, "1107": 5, "s2059798317007859": 5, "2014": 5, "compress": 5, "recarrai": 5, "enforc": 5, "regardless": 5, "matrix": 5, "voxel": 5, "voxel_s": 5, "offset": 5, "nxstart": 5, "nystart": 5, "nzstart": 5, "denot": 5, "map": 5, "unitcel": 5, "discret": 6, "parallel": 6, "essenti": 6, "intersect": 6, "anchor": 6, "resolut": 6, "algebra": 6, "wise": 6, "basi": 6, "consist": 6, "itself": 6, "equival": 6, "know": 6, "thu": 6, "simpli": 6, "subtract": 6, "two": 6, "get": 6, "multipli": 6, "combin": 6, "within": 6, "code": [6, 7, 8], "function": 6, "represent": 6, "abstract": 6, "straightforward": 6, "chimerax": 6, "disk": 6, "g55c0241": 7, "date": 7, "oct": 7, "30": 7, "2023": 7, "citat": 7, "easier": 7, "lesser": 7, "gnu": 7, "public": 7, "licens": 7, "distribut": 7, "github": 7, "repositori": 7, "mdanalysi": 7, "particip": 7, "ask": 7, "mdnalysi": 7, "discuss": 7, "mail": 7, "join": 7, "report": 7, "problem": 7, "enhanc": 7, "issu": 7, "tracker": 7, "open": 7, "welcom": 7, "fork": 7, "submit": 7, "manag": 8, "compil": 8, "eas": 8, "recommend": 8, "fulli": 8, "updat": 8, "environ": 8, "forg": 8, "channel": 8, "achiev": 8, "config": 8, "been": 8, "enabl": 8, "automat": 8, "download": 8, "appropri": 8, "upgrad": 8, "attempt": 8, "step": 8, "learn": 8, "switch": 8}, "objects": {"": [[6, 0, 0, "-", "gridData"]], "gridData": [[3, 0, 0, "-", "OpenDX"], [1, 0, 0, "-", "core"], [4, 0, 0, "-", "gOpenMol"], [5, 0, 0, "-", "mrc"]], "gridData.OpenDX": [[3, 1, 1, "", "DXInitObject"], [3, 3, 1, "", "DXParseError"], [3, 1, 1, "", "DXParser"], [3, 3, 1, "", "DXParserNoTokens"], [3, 1, 1, "", "DXclass"], [3, 1, 1, "", "array"], [3, 1, 1, "", "field"], [3, 1, 1, "", "gridconnections"], [3, 1, 1, "", "gridpositions"]], "gridData.OpenDX.DXInitObject": [[3, 2, 1, "", "initialize"]], "gridData.OpenDX.DXParser": [[3, 2, 1, "", "apply_parser"], [3, 2, 1, "", "parse"], [3, 2, 1, "", "set_parser"], [3, 2, 1, "", "use_parser"]], "gridData.OpenDX.DXclass": [[3, 2, 1, "", "ndformat"], [3, 2, 1, "", "write"]], "gridData.OpenDX.array": [[3, 4, 1, "", "dx_types"], [3, 4, 1, "", "np_types"], [3, 2, 1, "", "write"]], "gridData.OpenDX.field": [[3, 2, 1, "", "add"], [3, 2, 1, "", "add_comment"], [3, 2, 1, "", "histogramdd"], [3, 2, 1, "", "read"], [3, 2, 1, "", "sorted_components"], [3, 2, 1, "", "write"]], "gridData.OpenDX.gridconnections": [[3, 2, 1, "", "write"]], "gridData.OpenDX.gridpositions": [[3, 2, 1, "", "edges"], [3, 2, 1, "", "write"]], "gridData.core": [[1, 1, 1, "", "Grid"], [1, 6, 1, "", "ndmeshgrid"]], "gridData.core.Grid": [[1, 2, 1, "", "centers"], [1, 2, 1, "", "check_compatible"], [1, 4, 1, "", "default_format"], [1, 4, 1, "", "delta"], [1, 4, 1, "", "edges"], [1, 2, 1, "", "export"], [1, 4, 1, "", "grid"], [1, 5, 1, "", "interpolated"], [1, 5, 1, "", "interpolation_spline_order"], [1, 2, 1, "", "load"], [1, 4, 1, "", "metadata"], [1, 4, 1, "", "midpoints"], [1, 4, 1, "", "origin"], [1, 2, 1, "", "resample"], [1, 2, 1, "", "resample_factor"], [1, 2, 1, "", "save"]], "gridData.gOpenMol": [[4, 1, 1, "", "Plt"]], "gridData.gOpenMol.Plt": [[4, 5, 1, "", "edges"], [4, 2, 1, "", "histogramdd"], [4, 2, 1, "", "read"]], "gridData.mrc": [[5, 1, 1, "", "MRC"]], "gridData.mrc.MRC": [[5, 4, 1, "", "array"], [5, 4, 1, "", "delta"], [5, 5, 1, "", "edges"], [5, 4, 1, "", "header"], [5, 2, 1, "", "histogramdd"], [5, 4, 1, "", "origin"], [5, 4, 1, "", "rank"], [5, 2, 1, "", "read"], [5, 5, 1, "", "shape"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method", "3": "py:exception", "4": "py:attribute", "5": "py:property", "6": "py:function"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"], "3": ["py", "exception", "Python exception"], "4": ["py", "attribute", "Python attribute"], "5": ["py", "property", "Python property"], "6": ["py", "function", "Python function"]}, "titleterms": {"basic": 0, "us": 0, "load": 0, "data": [0, 4, 5, 6, 7], "write": [0, 3], "out": 0, "subtract": 0, "two": 0, "densiti": 0, "resampl": 0, "core": 1, "function": [1, 3], "store": 1, "n": 1, "d": 1, "grid": [1, 4, 6], "griddata": [1, 2, 6], "class": [1, 3, 4, 5], "format": [2, 4, 5, 6], "support": 2, "file": [2, 3, 4, 6], "avail": 2, "specif": 2, "modul": 2, "opendx": 3, "routin": 3, "read": [3, 6], "simpl": 3, "known": 3, "issu": 3, "build": 3, "dx": 3, "object": 3, "from": 3, "numpi": 3, "arrai": 3, "A": 3, "gopenmol": 4, "plt": 4, "background": 4, "binari": 4, "mrc": 5, "ccp4": 5, "volumetr": [5, 7], "handl": [6, 7], "overview": 6, "descript": 6, "construct": 6, "griddataformat": [7, 8], "python": 7, "instal": 8, "conda": 8, "pip": 8}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx.ext.viewcode": 1, "sphinx.ext.todo": 2, "sphinx": 60}, "alltitles": {"Basic use": [[0, "basic-use"]], "Loading data": [[0, "loading-data"]], "Writing out data": [[0, "writing-out-data"]], "Subtracting two densities": [[0, "subtracting-two-densities"]], "Resampling": [[0, "resampling"]], "Core functionality for storing n-D grids \u2014 gridData.core": [[1, "core-functionality-for-storing-n-d-grids-griddata-core"]], "Classes and functions": [[1, "classes-and-functions"], [3, "classes-and-functions"]], "Formats": [[2, "formats"], [6, "formats"]], "Supported file formats": [[2, "supported-file-formats"]], "Available file formats in gridData": [[2, "id2"]], "Format-specific modules": [[2, "format-specific-modules"]], "OpenDX \u2014 routines to read and write simple OpenDX files": [[3, "opendx-routines-to-read-and-write-simple-opendx-files"]], "Reading and writing OpenDX files": [[3, "reading-and-writing-opendx-files"]], "Known issues for writing OpenDX files": [[3, "known-issues-for-writing-opendx-files"]], "Building a dx object from a numpy array A": [[3, "building-a-dx-object-from-a-numpy-array-a"]], "Building a dx object from a dx file": [[3, "building-a-dx-object-from-a-dx-file"]], "gOpenMol \u2014 the gOpenMol plt format": [[4, "gopenmol-the-gopenmol-plt-format"]], "Background": [[4, "background"]], "Grid data plt file format": [[4, "grid-data-plt-file-format"]], "Format of the binary *.plt file": [[4, "format-of-the-binary-plt-file"]], "Binary *.plt (grid) file format": [[4, "binary-plt-grid-file-format"]], "Formatted *.plt (grid) file format": [[4, "formatted-plt-grid-file-format"]], "Classes": [[4, "classes"], [5, "classes"]], "mrc \u2014 the MRC/CCP4 volumetric data format": [[5, "mrc-the-mrc-ccp4-volumetric-data-format"]], "Handling grids of data \u2014 gridData": [[6, "handling-grids-of-data-griddata"]], "Overview": [[6, "overview"]], "Description": [[6, "description"]], "Reading grid data files": [[6, "reading-grid-data-files"]], "Constructing a Grid": [[6, "constructing-a-grid"]], "GridDataFormats: Handling volumetric data in Python": [[7, "griddataformats-handling-volumetric-data-in-python"]], "Installation": [[8, "installation"]], "Installing GridDataFormats with conda": [[8, "installing-griddataformats-with-conda"]], "Installing GridDataFormats with pip": [[8, "installing-griddataformats-with-pip"]]}, "indexentries": {"grid (class in griddata.core)": [[1, "gridData.core.Grid"]], "centers() (griddata.core.grid method)": [[1, "gridData.core.Grid.centers"]], "check_compatible() (griddata.core.grid method)": [[1, "gridData.core.Grid.check_compatible"]], "default_format (griddata.core.grid attribute)": [[1, "gridData.core.Grid.default_format"]], "delta (griddata.core.grid attribute)": [[1, "gridData.core.Grid.delta"]], "edges (griddata.core.grid attribute)": [[1, "gridData.core.Grid.edges"]], "export() (griddata.core.grid method)": [[1, "gridData.core.Grid.export"]], "grid (griddata.core.grid attribute)": [[1, "gridData.core.Grid.grid"]], "griddata.core": [[1, "module-gridData.core"]], "interpolated (griddata.core.grid property)": [[1, "gridData.core.Grid.interpolated"]], "interpolation_spline_order (griddata.core.grid property)": [[1, "gridData.core.Grid.interpolation_spline_order"]], "load() (griddata.core.grid method)": [[1, "gridData.core.Grid.load"]], "metadata (griddata.core.grid attribute)": [[1, "gridData.core.Grid.metadata"]], "midpoints (griddata.core.grid attribute)": [[1, "gridData.core.Grid.midpoints"]], "module": [[1, "module-gridData.core"], [3, "module-gridData.OpenDX"], [4, "module-gridData.gOpenMol"], [5, "module-gridData.mrc"], [6, "module-gridData"]], "ndmeshgrid() (in module griddata.core)": [[1, "gridData.core.ndmeshgrid"]], "origin (griddata.core.grid attribute)": [[1, "gridData.core.Grid.origin"]], "resample() (griddata.core.grid method)": [[1, "gridData.core.Grid.resample"]], "resample_factor() (griddata.core.grid method)": [[1, "gridData.core.Grid.resample_factor"]], "save() (griddata.core.grid method)": [[1, "gridData.core.Grid.save"]], "dxinitobject (class in griddata.opendx)": [[3, "gridData.OpenDX.DXInitObject"]], "dxparseerror": [[3, "gridData.OpenDX.DXParseError"]], "dxparser (class in griddata.opendx)": [[3, "gridData.OpenDX.DXParser"]], "dxparsernotokens": [[3, "gridData.OpenDX.DXParserNoTokens"]], "dxclass (class in griddata.opendx)": [[3, "gridData.OpenDX.DXclass"]], "add() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.add"]], "add_comment() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.add_comment"]], "apply_parser() (griddata.opendx.dxparser method)": [[3, "gridData.OpenDX.DXParser.apply_parser"]], "array (class in griddata.opendx)": [[3, "gridData.OpenDX.array"]], "dx_types (griddata.opendx.array attribute)": [[3, "gridData.OpenDX.array.dx_types"]], "edges() (griddata.opendx.gridpositions method)": [[3, "gridData.OpenDX.gridpositions.edges"]], "field (class in griddata.opendx)": [[3, "gridData.OpenDX.field"]], "griddata.opendx": [[3, "module-gridData.OpenDX"]], "gridconnections (class in griddata.opendx)": [[3, "gridData.OpenDX.gridconnections"]], "gridpositions (class in griddata.opendx)": [[3, "gridData.OpenDX.gridpositions"]], "histogramdd() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.histogramdd"]], "initialize() (griddata.opendx.dxinitobject method)": [[3, "gridData.OpenDX.DXInitObject.initialize"]], "ndformat() (griddata.opendx.dxclass method)": [[3, "gridData.OpenDX.DXclass.ndformat"]], "np_types (griddata.opendx.array attribute)": [[3, "gridData.OpenDX.array.np_types"]], "parse() (griddata.opendx.dxparser method)": [[3, "gridData.OpenDX.DXParser.parse"]], "read() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.read"]], "set_parser() (griddata.opendx.dxparser method)": [[3, "gridData.OpenDX.DXParser.set_parser"]], "sorted_components() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.sorted_components"]], "use_parser() (griddata.opendx.dxparser method)": [[3, "gridData.OpenDX.DXParser.use_parser"]], "write() (griddata.opendx.dxclass method)": [[3, "gridData.OpenDX.DXclass.write"]], "write() (griddata.opendx.array method)": [[3, "gridData.OpenDX.array.write"]], "write() (griddata.opendx.field method)": [[3, "gridData.OpenDX.field.write"]], "write() (griddata.opendx.gridconnections method)": [[3, "gridData.OpenDX.gridconnections.write"]], "write() (griddata.opendx.gridpositions method)": [[3, "gridData.OpenDX.gridpositions.write"]], "plt (class in griddata.gopenmol)": [[4, "gridData.gOpenMol.Plt"]], "edges (griddata.gopenmol.plt property)": [[4, "gridData.gOpenMol.Plt.edges"]], "griddata.gopenmol": [[4, "module-gridData.gOpenMol"]], "histogramdd() (griddata.gopenmol.plt method)": [[4, "gridData.gOpenMol.Plt.histogramdd"]], "read() (griddata.gopenmol.plt method)": [[4, "gridData.gOpenMol.Plt.read"]], "mrc (class in griddata.mrc)": [[5, "gridData.mrc.MRC"]], "array (griddata.mrc.mrc attribute)": [[5, "gridData.mrc.MRC.array"]], "delta (griddata.mrc.mrc attribute)": [[5, "gridData.mrc.MRC.delta"]], "edges (griddata.mrc.mrc property)": [[5, "gridData.mrc.MRC.edges"]], "griddata.mrc": [[5, "module-gridData.mrc"]], "header (griddata.mrc.mrc attribute)": [[5, "gridData.mrc.MRC.header"]], "histogramdd() (griddata.mrc.mrc method)": [[5, "gridData.mrc.MRC.histogramdd"]], "origin (griddata.mrc.mrc attribute)": [[5, "gridData.mrc.MRC.origin"]], "rank (griddata.mrc.mrc attribute)": [[5, "gridData.mrc.MRC.rank"]], "read() (griddata.mrc.mrc method)": [[5, "gridData.mrc.MRC.read"]], "shape (griddata.mrc.mrc property)": [[5, "gridData.mrc.MRC.shape"]], "griddata": [[6, "module-gridData"]]}})
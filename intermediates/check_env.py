import sys
import importlib
print("Python:", sys.version)
for lib in ["pandas", "numpy", "scipy", "sklearn", "anndata", "scanpy", "umap"]:
    try:
        m = importlib.import_module(lib)
        v = getattr(m, "__version__", "?")
        print(f"  {lib}: {v}")
    except ImportError:
        print(f"  {lib}: NOT AVAILABLE")

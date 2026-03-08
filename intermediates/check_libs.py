import sys
print("Python:", sys.version)
libs = ["pandas","numpy","scipy","sklearn","umap","anndata","scanpy",
        "statsmodels","matplotlib","seaborn","networkx","igraph"]
for lib in libs:
    try:
        import importlib
        m = importlib.import_module(lib)
        print(f"  OK   {lib}: {getattr(m,'__version__','?')}")
    except ImportError as e:
        print(f"  MISS {lib}: {e}")

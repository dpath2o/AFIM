# tools/find_nested_defs.py
import ast, sys, pathlib

def find_nested_defs(py_path: str):
    src = pathlib.Path(py_path).read_text()
    tree = ast.parse(src)
    parents = []

    class V(ast.NodeVisitor):
        def visit_FunctionDef(self, node: ast.FunctionDef):
            parents.append(node)
            self.generic_visit(node)
            parents.pop()

        def visit_AsyncFunctionDef(self, node):
            self.visit_FunctionDef(node)

        def visit(self, node):
            # report any function nested within another function
            if isinstance(node, ast.FunctionDef) and parents and parents[-1] is not node:
                parent = parents[-1]
                print(f"{py_path}:{node.lineno} nested def '{node.name}' inside '{parent.name}'")
            super().visit(node)

    V().visit(tree)

if __name__ == "__main__":
    for p in sys.argv[1:] or ["src/sea_ice_toolbox.py"]:
        find_nested_defs(p)

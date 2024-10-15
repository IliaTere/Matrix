import difflib
import re
import clang.cindex

clang.cindex.Config.set_library_file('/opt/homebrew/Cellar/llvm/19.1.1/lib/libclang.dylib')

def read_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.read()

def compare_text(text1, text2):
    matcher = difflib.SequenceMatcher(None, text1, text2)
    similarity = matcher.ratio()
    return similarity

def compare_ast(code1, code2):
    index = clang.cindex.Index.create()
    tu1 = index.parse('temp1.cpp', args=['-std=c++11'], unsaved_files=[('temp1.cpp', code1)])
    tu2 = index.parse('temp2.cpp', args=['-std=c++11'], unsaved_files=[('temp2.cpp', code2)])
    return ast_similarity(tu1.cursor, tu2.cursor)

def ast_similarity(cursor1, cursor2):
    if cursor1.kind != cursor2.kind:
        return False
    if cursor1.spelling != cursor2.spelling:
        return False
    if len(list(cursor1.get_children())) != len(list(cursor2.get_children())):
        return False
    for child1, child2 in zip(cursor1.get_children(), cursor2.get_children()):
        if not ast_similarity(child1, child2):
            return False
    return True

def compare_identifiers(code1, code2):
    identifiers1 = set(re.findall(r'\b\w+\b', code1))
    identifiers2 = set(re.findall(r'\b\w+\b', code2))
    common_identifiers = identifiers1 & identifiers2
    return len(common_identifiers) / max(len(identifiers1), len(identifiers2))

def compare_comments(code1, code2):
    comments1 = re.findall(r'//.*|/\*.*?\*/', code1, re.DOTALL)
    comments2 = re.findall(r'//.*|/\*.*?\*/', code2, re.DOTALL)
    common_comments = set(comments1) & set(comments2)
    return len(common_comments) / max(len(comments1), len(comments2))

def compare_style(code1, code2):
    style_pattern = re.compile(r'\s+')
    style1 = style_pattern.findall(code1)
    style2 = style_pattern.findall(code2)
    return style1 == style2

def compare_classes_and_functions(code1, code2):
    classes_and_functions1 = set(re.findall(r'\b(class|struct|void|int|float|double|bool|char|std::string|std::vector)\s+\w+\b', code1))
    classes_and_functions2 = set(re.findall(r'\b(class|struct|void|int|float|double|bool|char|std::string|std::vector)\s+\w+\b', code2))
    common_elements = classes_and_functions1 & classes_and_functions2
    return len(common_elements) / max(len(classes_and_functions1), len(classes_and_functions2))

def compare_namespaces(code1, code2):
    namespaces1 = set(re.findall(r'\bnamespace\s+\w+\b', code1))
    namespaces2 = set(re.findall(r'\bnamespace\s+\w+\b', code2))
    common_namespaces = namespaces1 & namespaces2
    return len(common_namespaces) / max(len(namespaces1), len(namespaces2))

def compare_macros_and_constants(code1, code2):
    macros_and_constants1 = set(re.findall(r'\b(const|#define)\s+\w+\b', code1))
    macros_and_constants2 = set(re.findall(r'\b(const|#define)\s+\w+\b', code2))
    
    if not macros_and_constants1 and not macros_and_constants2:
        return 0.0
    
    common_elements = macros_and_constants1 & macros_and_constants2
    return len(common_elements) / max(len(macros_and_constants1), len(macros_and_constants2))
def compare_libraries_and_headers(code1, code2):
    libraries_and_headers1 = set(re.findall(r'#include\s+[<"]\w+\.h[>"]', code1))
    libraries_and_headers2 = set(re.findall(r'#include\s+[<"]\w+\.h[>"]', code2))
    common_elements = libraries_and_headers1 & libraries_and_headers2
    return len(common_elements) / max(len(libraries_and_headers1), len(libraries_and_headers2))

def main(file1, file2):
    text1 = read_file(file1)
    text2 = read_file(file2)

    text_similarity = compare_text(text1, text2)
    ast_similarity = compare_ast(text1, text2)
    identifier_similarity = compare_identifiers(text1, text2)
    comment_similarity = compare_comments(text1, text2)
    style_similarity = compare_style(text1, text2)
    class_function_similarity = compare_classes_and_functions(text1, text2)
    namespace_similarity = compare_namespaces(text1, text2)
    macro_constant_similarity = compare_macros_and_constants(text1, text2)
    library_header_similarity = compare_libraries_and_headers(text1, text2)

    print(f"Text Similarity: {text_similarity:.2f}")
    print(f"AST Similarity: {ast_similarity}")
    print(f"Identifier Similarity: {identifier_similarity:.2f}")
    print(f"Comment Similarity: {comment_similarity:.2f}")
    print(f"Style Similarity: {style_similarity}")
    print(f"Class and Function Similarity: {class_function_similarity:.2f}")
    print(f"Namespace Similarity: {namespace_similarity:.2f}")
    print(f"Macro and Constant Similarity: {macro_constant_similarity:.2f}")
    print(f"Library and Header Similarity: {library_header_similarity:.2f}")

if __name__ == "__main__":
    file1 = "qwe.cpp"
    file2 = "solve.hpp"
    main(file1, file2)
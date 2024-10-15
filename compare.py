import difflib
import re
import clang.cindex
clang.cindex.Config.set_library_file('/usr/lib/x86_64-linux-gnu/libclang-14.so')
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

def main(file1, file2):
    text1 = read_file(file1)
    text2 = read_file(file2)

    text_similarity = compare_text(text1, text2)
    ast_similarity = compare_ast(text1, text2)
    identifier_similarity = compare_identifiers(text1, text2)
    comment_similarity = compare_comments(text1, text2)
    style_similarity = compare_style(text1, text2)

    print(f"Text Similarity: {text_similarity:.2f}")
    print(f"AST Similarity: {ast_similarity}")
    print(f"Identifier Similarity: {identifier_similarity:.2f}")
    print(f"Comment Similarity: {comment_similarity:.2f}")
    print(f"Style Similarity: {style_similarity}")

if __name__ == "__main__":
    file1 = "qwe.cpp"
    file2 = "solve.hpp"
    main(file1, file2)
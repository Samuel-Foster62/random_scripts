#!/usr/bin/env python3

import sys
from bs4 import BeautifulSoup

def parse_tables(html_file):
    with open(html_file, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")

    tables = soup.find_all("table")
    parsed = []

    for t_idx, table in enumerate(tables):
        rows = table.find_all("tr")
        table_data = []

        for row in rows:
            cols = row.find_all(["td", "th"])
            text = [c.get_text(strip=True) for c in cols]
            if text:
                table_data.append(text)

        if table_data:
            parsed.append(table_data)

    return parsed


def preview_tables(tables):
    print("\nDetected tables:\n")

    for i, table in enumerate(tables):
        print(f"Table {i}:")
        for row in table[:5]:  # first 5 rows
            print("  " + " | ".join(row))
        print("-" * 60)


def find_image_template_column(headers):
    for i, h in enumerate(headers):
        if "image template" in h.lower():
            return i
    return None


def extract_image_templates(table):
    headers = table[0]
    col_idx = find_image_template_column(headers)

    if col_idx is None:
        print("No 'Image template' column found in this table.")
        return

    print("\nImage templates:\n")
    for row in table[1:]:
        if len(row) > col_idx:
            print(row[col_idx])


def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_templates.py <html_file>")
        sys.exit(1)

    html_file = sys.argv[1]

    tables = parse_tables(html_file)

    if not tables:
        print("No tables found.")
        return

    preview_tables(tables)

    # user selection
    while True:
        try:
            choice = int(input("\nSelect table number: "))
            if 0 <= choice < len(tables):
                break
            else:
                print("Invalid number.")
        except ValueError:
            print("Please enter a valid integer.")

    extract_image_templates(tables[choice])


if __name__ == "__main__":
    main()
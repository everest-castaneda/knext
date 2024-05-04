from collections import defaultdict


def _parse_entries(root):
    entry_dict = defaultdict(list)
    for entries in root.findall('entry'):
        for key, items in entries.attrib.items():
            entry_dict[key].append(items)

    entry_id=[]
    entry_name=[]
    entry_type=[]
    for key, items in entry_dict.items():
        if key == 'id':
            for i in items:
                entry_id.append(i)
        if key == 'name':
            for i in items:
                entry_name.append(i)
        if key == 'type':
            for i in items:
                entry_type.append(i)

    return entry_id, entry_name, entry_type
import jsonschema
import json
import os
import sys

with open('test_schema.json', 'r') as f:
    schema_data = f.read()
schema = json.loads(schema_data)

with open("test.json", "r") as f:
    json_data = f.read()
tests = json.loads(json_data)

for test_name, test_value in tests.items():
    # try:
    jsonschema.validate(test_value, schema)
    # except:
    #     print("Error: test {} is incorrectly specified in test.json".format(
    #         key))
    #     raise SystemExit

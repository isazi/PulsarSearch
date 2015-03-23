#!/usr/bin/env python
# Copyright 2015 Alessio Sclocco <a.sclocco@vu.nl>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

def get_tables(queue):
    """Get a list of the tables"""
    queue.execute("SHOW TABLES")
    return queue.fetchall()

def create_table(queue, table):
    """Create a table to store the pipeline performance measurements."""
    queue.execute("CREATE table " + table + "(id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT, DMs INTEGER NOT NULL, periods INTEGER NOT NULL, searchTime FLOAT UNSIGNED NOT NULL, inputHandling FLOAT UNSIGNED NOT NULL, dedispersion FLOAT UNSIGNED NOT NULL, transpose FLOAT UNSIGNED NOT NULL, snrD FLOAT UNSIGNED NOT NULL, folding FLOAT UNSIGNED NOT NULL, snrF FLOAT UNSIGNED NOT NULL, outputCopy FLOAT UNSIGNED NOT NULL)")

def delete_table(queue, table):
    """Delete table."""
    queue.execute("DROP table " + table)

def load_file(queue, table, input_file):
    """Load input_file into a table in the database."""
    for line in input_file:
        if (line[0] != "#") and (line[0] != "\n"):
            items = line.split(sep=" ")
            queue.execute("INSERT INTO " + table + " VALUES (NULL, " + items[0] + ", " + items[1] + ", " + items[2] + ", " + items[3] + ", " + items[6] + ", " + items[9] + ", " + items[12] + ", " + items[15] + ", " + items[18] + ", " + str(float(items[21]) + float(items[24])) + ")")

def print_results(confs):
    """Print the result tuples."""
    for conf in confs:
        for internalConf in conf:
            for item in internalConf:
                print(item, end=" ")
            print()
        print("\n")

def get_dm_range(queue, table):
    """Return the DMs used in the scenario."""
    queue.execute("SELECT DISTINCT DMs FROM " + table + " ORDER BY DMs")
    return queue.fetchall()

def get_period_range(queue, table, dm):
    """Return the periods used in the scenario."""
    queue.execute("SELECT DISTINCT periods FROM " + table + " WHERE (DMs = " + dm + ") ORDER BY periods")
    return queue.fetchall()


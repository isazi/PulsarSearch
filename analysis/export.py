#!/usr/bin/env python
# Copyright 2014 Alessio Sclocco <a.sclocco@vu.nl>
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

import manage

def print_time(queue, table):
    confs = list()
    dms_range = manage.get_dm_range(queue, table)
    for dm in dms_range:
        internalConfs = list()
        period_range = manage.get_period_range(queue, table, str(dm[0]))
        for period in period_range:
            queue.execute("SELECT searchTime,inputHandling,dedispersion,transpose,snrD,folding,snrF,outputCopy FROM " + table + " WHERE (DMs = " + str(dm[0]) + " AND periods = " + str(period[0]) + ")")
            best = queue.fetchall()
            internalConfs.append([dm[0], period[0], best[0][0], best[0][1], best[0][2], best[0][3], best[0][4], best[0][5], best[0][6], best[0][7]])
        confs.append(internalConfs)
    return confs


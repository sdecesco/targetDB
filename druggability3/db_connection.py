#!/usr/bin/env python
import MySQLdb
import itertools


def open_db(database,pwd=None,user=None):

    db = MySQLdb.connect("localhost", user, pwd, database,charset='utf8',use_unicode=True)
    db.query("""SET GROUP_CONCAT_MAX_LEN = 30000""")
    db_obj = db_connector(db,database,pwd,user)
    return db_obj


def close_db(db_obj):
    try:
        db_obj.db.close()
        return True
    except:
        return False


class db_connector:
    def __init__(self, db,database,pwd,user):
        self.db = db
        self.database = database
        self.pwd = pwd
        self.user = user

    def connect(self,database):
        self.db=MySQLdb.connect("localhost", self.user, self.pwd, database,charset='utf8',use_unicode=True)
        self.db.query("""SET GROUP_CONCAT_MAX_LEN = 30000""")

    def get(self, query):
        cursor_collect = self.db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        try:
            cursor_collect.execute(query)
        except (AttributeError,MySQLdb.OperationalError):
            self.connect(self.database)
            cursor_collect = self.db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
            cursor_collect.execute(query)
        results = cursor_collect.fetchall()
        cursor_collect.close()
        return results

    def write(self, query, values):
        cursor_write = self.db.cursor()
        cursor_write.execute(query, values)
        self.db.commit()
        cursor_write.close()
        return True

    def multi_write(self, query, values, update=''):
        if not values:
            return False
        cursor_write = self.db.cursor()
        n_column = len(query.split('(')[1].split(','))
        tup_value = tuple(itertools.chain(*values))
        query_2 = "VALUES"
        values_base = "("
        for i in range(0, n_column):
            values_base += "%s,"
        values_base = values_base.rstrip(',')
        values_base += '),'
        if update:
            update = ' ' + update
        count = 0
        limit = 500
        if len(values) <= limit:
            limit = len(values)
        n_loop = 0
        excess = len(values)%limit
        for i in range(len(values)):
            query_2 += values_base
            count += 1
            if count % limit == 0:
                values_to_query = tup_value[(n_loop * limit) * n_column:((n_loop * limit) * n_column) + (limit * n_column)]
                query_2 = query_2.rstrip(',')
                full_query = query + query_2 + update
                cursor_write.execute(full_query, values_to_query)
                self.db.commit()
                query_2 = "VALUES"
                count = 0
                n_loop += 1
        if excess != 0:
            values_to_query = tup_value[(n_loop * limit) * n_column:((n_loop * limit) * n_column) + (excess * n_column)]
            query_2 = query_2.rstrip(',')
            full_query = query + query_2 + update
            cursor_write.execute(full_query, values_to_query)
            self.db.commit()


        cursor_write.close()
        return True

    def close(self):
        try:
            self.db.close()
            return True
        except:
            return False


if __name__ == "__main__":
    pass

import time


def time_of_function(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        execution_time = round(time.time() - start_time, 3)
        print(f'Время выполнения: {execution_time} сек.')
        return result
    return wrapper

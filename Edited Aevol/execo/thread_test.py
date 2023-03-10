from threading import Thread

def thread_main():
    print "hello"
    sleep 100000

if __name__ == "__main__":
    threads = {}
    for i in range(0,100):
        t = Thread(target=thread_main)
        threads[t] = {'host': i}
        t.daemon = True
        t.start()
        sleep(3)

import time as t

class MyTimer():

    def __init__(self):
        self.unit=['year','month','day','H','M','S']
        self.prompt="not begin yet"
        self.lasted=[]
        self.begin=0
        self.end=0


    #开始计时
    def __str__(self):
        return self.prompt
    __repr__=__str__

    def __add__(self,other):
        prompt="run time:"
        result=[]
        for index in range(6):
            result.append(self.lasted[index]+other.lasted[index])
            if result[index]:
                prompt+=(str(result[index])+self.unit[index])
        return prompt

    
    def start(self,string=None):
        out = "..."
        if string != None:
            out = string+","+out
        self.begin=t.localtime()
        self.prompt="please call stop() to shut down the timer"
        print(out)

    #停止计时
    def stop(self):
        if not self.begin:
            print("please call start() to start up the timer")
        else:
          self.end=t.localtime()
          self._calc()
          print(str(self.prompt))


    #内部方法
    def _calc(self):
        self.lasted=[]
        self.prompt="run time:"
        for index in range(6):
            self.lasted.append(self.end[index]-self.begin[index])
            if self.lasted[index]:
               self.prompt+=(str(self.lasted[index])+self.unit[index])
        #为下一轮计时初始化变量
        self.begin=0
        self.end=0

    
        

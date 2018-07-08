c    ------------------
      function cnm(nn)
c    ------------------
      character(len=3) cnm
      character ch1,ch2,ch3

      ch1 = char(int(mod(nn,1000)/100)+48)
      ch2 = char(int(mod(nn,100)/10)+48)
      ch3 = char(int(mod(nn,10))+48)
      cnm = ch1//ch2//ch3

      return
      end function cnm

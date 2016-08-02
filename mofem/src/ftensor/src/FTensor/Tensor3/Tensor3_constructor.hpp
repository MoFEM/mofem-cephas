/* A helper class that allows simple initialization of the Tensor3,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
class Tensor3_constructor;

template<class T>
class Tensor3_constructor<T,2,2,2>
{
public:
  Tensor3_constructor(T data[2][2][2], T d000, T d001, T d010, T d011,
			 T d100, T d101, T d110, T d111)
  {
    data[0][0][0]=d000;
    data[0][0][1]=d001;
    data[0][1][0]=d010;
    data[0][1][1]=d011;
    data[1][0][0]=d100;
    data[1][0][1]=d101;
    data[1][1][0]=d110;
    data[1][1][1]=d111;
  }
};

template<class T>
class Tensor3_constructor<T,3,3,3>
{
public:
  Tensor3_constructor(T data[3][3][3], T d000, T d001, T d002,
                      T d010, T d011, T d012,
                      T d020, T d021, T d022,
                      T d100, T d101, T d102,
                      T d110, T d111, T d112,
                      T d120, T d121, T d122,
                      T d200, T d201, T d202,
                      T d210, T d211, T d212,
                      T d220, T d221, T d222)
  {
    data[0][0][0]=d000;
    data[0][0][1]=d001;
    data[0][0][2]=d002;

    data[0][1][0]=d010;
    data[0][1][1]=d011;
    data[0][1][2]=d012;

    data[0][2][0]=d020;
    data[0][2][1]=d021;
    data[0][2][2]=d022;

    data[1][0][0]=d100;
    data[1][0][1]=d101;
    data[1][0][2]=d102;

    data[1][1][0]=d110;
    data[1][1][1]=d111;
    data[1][1][2]=d112;

    data[1][2][0]=d120;
    data[1][2][1]=d121;
    data[1][2][2]=d122;

    data[2][0][0]=d200;
    data[2][0][1]=d201;
    data[2][0][2]=d202;

    data[2][1][0]=d210;
    data[2][1][1]=d211;
    data[2][1][2]=d212;

    data[2][2][0]=d220;
    data[2][2][1]=d221;
    data[2][2][2]=d222;
  }
};

template<class T>
class Tensor3_constructor<T,4,4,4>
{
public:
  Tensor3_constructor(T data[4][4][4], T d000, T d001, T d002, T d003,
                      T d010, T d011, T d012, T d013,
                      T d020, T d021, T d022, T d023,
                      T d030, T d031, T d032, T d033,

                      T d100, T d101, T d102, T d103,
                      T d110, T d111, T d112, T d113,
                      T d120, T d121, T d122, T d123,
                      T d130, T d131, T d132, T d133,

                      T d200, T d201, T d202, T d203,
                      T d210, T d211, T d212, T d213,
                      T d220, T d221, T d222, T d223,
                      T d230, T d231, T d232, T d233,

                      T d300, T d301, T d302, T d303,
                      T d310, T d311, T d312, T d313,
                      T d320, T d321, T d322, T d323,
                      T d330, T d331, T d332, T d333)
  {
    data[0][0][0]=d000;
    data[0][0][1]=d001;
    data[0][0][2]=d002;
    data[0][0][3]=d003;

    data[0][1][0]=d010;
    data[0][1][1]=d011;
    data[0][1][2]=d012;
    data[0][1][3]=d013;

    data[0][2][0]=d020;
    data[0][2][1]=d021;
    data[0][2][2]=d022;
    data[0][2][3]=d023;

    data[0][3][0]=d030;
    data[0][3][1]=d031;
    data[0][3][2]=d032;
    data[0][3][3]=d033;

    data[1][0][0]=d100;
    data[1][0][1]=d101;
    data[1][0][2]=d102;
    data[1][0][3]=d103;

    data[1][1][0]=d110;
    data[1][1][1]=d111;
    data[1][1][2]=d112;
    data[1][1][3]=d113;

    data[1][2][0]=d120;
    data[1][2][1]=d121;
    data[1][2][2]=d122;
    data[1][2][3]=d123;

    data[1][3][0]=d130;
    data[1][3][1]=d131;
    data[1][3][2]=d132;
    data[1][3][3]=d133;

    data[2][0][0]=d200;
    data[2][0][1]=d201;
    data[2][0][2]=d202;
    data[2][0][3]=d203;

    data[2][1][0]=d210;
    data[2][1][1]=d211;
    data[2][1][2]=d212;
    data[2][1][3]=d213;

    data[2][2][0]=d220;
    data[2][2][1]=d221;
    data[2][2][2]=d222;
    data[2][2][3]=d223;

    data[2][3][0]=d230;
    data[2][3][1]=d231;
    data[2][3][2]=d232;
    data[2][3][3]=d233;

    data[3][0][0]=d300;
    data[3][0][1]=d301;
    data[3][0][2]=d302;
    data[3][0][3]=d303;

    data[3][1][0]=d310;
    data[3][1][1]=d311;
    data[3][1][2]=d312;
    data[3][1][3]=d313;

    data[3][2][0]=d320;
    data[3][2][1]=d321;
    data[3][2][2]=d322;
    data[3][2][3]=d323;

    data[3][3][0]=d330;
    data[3][3][1]=d331;
    data[3][3][2]=d332;
    data[3][3][3]=d333;
  }
};

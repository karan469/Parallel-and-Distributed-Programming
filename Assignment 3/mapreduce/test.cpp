#include <bits/stdc++.h>
using namespace std;
#define print(x) cout<<x<<'\n'

namespace first{
int a = 6;
void func(){
	print("Hello world!");
}

}

int main(int argc, char const *argv[])
{
	first::func();
	print(first::a);
	return 0;
}
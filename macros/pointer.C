void pointer()
{
  Double_t ppt[5] = { 1., 2., 3., 4., 5.};

  Double_t *pt[5];

  for (UInt_t i = 0; i < 5; i++ ) {

    pt[i] = &ppt[i];

    cout << "  ppt[" << i << "]" << ppt[i]
	 << " &ppt[" << i << "]" << &ppt[i]
	 << "   pt[" << i << "]" << pt[i]
	 << "   pt[" << i << "]" << *pt[i]
	 << endl;

    cout << " after --> *10" << endl;

    ppt[i] *= 10.;

    cout << "  ppt[" << i << "]" << ppt[i]
	 << " &ppt[" << i << "]" << &ppt[i]
	 << "   pt[" << i << "]" << pt[i]
	 << "   pt[" << i << "]" << *pt[i]
	 << endl;

    cout << " ------------------------" << endl;
  }

 //  ppt[0]1 &ppt[0]0x7fff2529dfb0   pt[0]0x7fff2529dfb0   pt[0]1
 // after --> *10
 //  ppt[0]10 &ppt[0]0x7fff2529dfb0   pt[0]0x7fff2529dfb0   pt[0]10
 // ------------------------
 //  ppt[1]2 &ppt[1]0x7fff2529dfb8   pt[1]0x7fff2529dfb8   pt[1]2
 // after --> *10
 //  ppt[1]20 &ppt[1]0x7fff2529dfb8   pt[1]0x7fff2529dfb8   pt[1]20
 // ------------------------
 //  ppt[2]3 &ppt[2]0x7fff2529dfc0   pt[2]0x7fff2529dfc0   pt[2]3
 // after --> *10
 //  ppt[2]30 &ppt[2]0x7fff2529dfc0   pt[2]0x7fff2529dfc0   pt[2]30
 // ------------------------
 //  ppt[3]4 &ppt[3]0x7fff2529dfc8   pt[3]0x7fff2529dfc8   pt[3]4
 // after --> *10
 //  ppt[3]40 &ppt[3]0x7fff2529dfc8   pt[3]0x7fff2529dfc8   pt[3]40
 // ------------------------
 //  ppt[4]5 &ppt[4]0x7fff2529dfd0   pt[4]0x7fff2529dfd0   pt[4]5
 // after --> *10
 //  ppt[4]50 &ppt[4]0x7fff2529dfd0   pt[4]0x7fff2529dfd0   pt[4]50
 // ------------------------

}

void pointer0()
{
  Double_t ppt[5] = { 1., 2., 3., 4., 5.};

  Double_t *pt = new Double_t[5];

  for (UInt_t i = 0; i < 5; i++ ) {

    pt[i] = ppt[i];
    //pt[i] = ppt[i];

    cout << "  ppt[" << i << "]" << ppt[i]
	 << " &ppt[" << i << "]" << &ppt[i]
	 << "   pt[" << i << "]" << pt[i]
	 << "  &pt[" << i << "]" << &pt[i]
	 << endl;

    cout << " after --> *10" << endl;

    ppt[i] *= 10.;

    cout << "  ppt[" << i << "]" << ppt[i]
	 << " &ppt[" << i << "]" << &ppt[i]
	 << "   pt[" << i << "]" << pt[i]
	 << "  &pt[" << i << "]" << &pt[i]
	 << endl;

    cout << " ------------------------" << endl;
  }

 //  ppt[0]1 &ppt[0]0x7fff2529eab0   pt[0]1  &pt[0]0x2112770
 // after --> *10
 //  ppt[0]10 &ppt[0]0x7fff2529eab0   pt[0]1  &pt[0]0x2112770
 // ------------------------
 //  ppt[1]2 &ppt[1]0x7fff2529eab8   pt[1]2  &pt[1]0x2112778
 // after --> *10
 //  ppt[1]20 &ppt[1]0x7fff2529eab8   pt[1]2  &pt[1]0x2112778
 // ------------------------
 //  ppt[2]3 &ppt[2]0x7fff2529eac0   pt[2]3  &pt[2]0x2112780
 // after --> *10
 //  ppt[2]30 &ppt[2]0x7fff2529eac0   pt[2]3  &pt[2]0x2112780
 // ------------------------
 //  ppt[3]4 &ppt[3]0x7fff2529eac8   pt[3]4  &pt[3]0x2112788
 // after --> *10
 //  ppt[3]40 &ppt[3]0x7fff2529eac8   pt[3]4  &pt[3]0x2112788
 // ------------------------
 //  ppt[4]5 &ppt[4]0x7fff2529ead0   pt[4]5  &pt[4]0x2112790
 // after --> *10
 //  ppt[4]50 &ppt[4]0x7fff2529ead0   pt[4]5  &pt[4]0x2112790
 // ------------------------


}

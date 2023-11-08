package com.hartwig.hmftools.errorprofile.recal;

import com.hartwig.hmftools.errorprofile.TandemRepeat;
import com.hartwig.hmftools.errorprofile.TandemRepeatFinder;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;

import org.junit.Test;

import junit.framework.TestCase;

public class TileBaseQualRecalTest
{
    @Test
    public void testComputeAdjustment()
    {
        String flowcell = "A00624:8:HHKYHDSXX:1";
        short lane = 1;
        short tile = 1182;
        boolean firstOfPair = true;
        byte rawBaseQuality = 37;
        TileAdjustmentKey tileAdjustmentKey = new TileAdjustmentKey(flowcell, lane, tile, firstOfPair, rawBaseQuality);
        TileBaseQualStats tileBaseQualStats = new TileBaseQualStats(tileAdjustmentKey);
        //tileBaseQualStats.addToCount(1, );

        // add some raw data
        tileBaseQualStats.addToCount(0, 10, 14936);
        tileBaseQualStats.addToCount(1, 1, 15120);
        tileBaseQualStats.addToCount(2, 4, 15177);
        tileBaseQualStats.addToCount(3, 5, 15186);
        tileBaseQualStats.addToCount(4, 2, 15208);
        tileBaseQualStats.addToCount(5, 4, 15229);
        tileBaseQualStats.addToCount(6, 3, 15237);
        tileBaseQualStats.addToCount(7, 3, 15266);
        tileBaseQualStats.addToCount(8, 1, 15248);
        tileBaseQualStats.addToCount(9, 0, 15224);
        tileBaseQualStats.addToCount(10, 2, 15252);
        tileBaseQualStats.addToCount(11, 3, 15239);
        tileBaseQualStats.addToCount(12, 1, 15257);
        tileBaseQualStats.addToCount(13, 2, 15217);
        tileBaseQualStats.addToCount(14, 1, 15244);
        tileBaseQualStats.addToCount(15, 2, 15260);
        tileBaseQualStats.addToCount(16, 0, 15196);
        tileBaseQualStats.addToCount(17, 1, 15197);
        tileBaseQualStats.addToCount(18, 2, 15230);
        tileBaseQualStats.addToCount(19, 1, 15211);
        tileBaseQualStats.addToCount(20, 1, 15224);
        tileBaseQualStats.addToCount(21, 2, 15185);
        tileBaseQualStats.addToCount(22, 1, 15175);
        tileBaseQualStats.addToCount(23, 0, 15206);
        tileBaseQualStats.addToCount(24, 2, 15183);
        tileBaseQualStats.addToCount(25, 0, 15179);
        tileBaseQualStats.addToCount(26, 1, 15147);
        tileBaseQualStats.addToCount(27, 1, 15169);
        tileBaseQualStats.addToCount(28, 3, 15143);
        tileBaseQualStats.addToCount(29, 1, 15169);
        tileBaseQualStats.addToCount(30, 1, 15150);
        tileBaseQualStats.addToCount(31, 1, 15121);
        tileBaseQualStats.addToCount(32, 1, 15143);
        tileBaseQualStats.addToCount(33, 1, 15155);
        tileBaseQualStats.addToCount(34, 1, 15106);
        tileBaseQualStats.addToCount(35, 2, 15139);
        tileBaseQualStats.addToCount(36, 1, 15092);
        tileBaseQualStats.addToCount(37, 1, 15081);
        tileBaseQualStats.addToCount(38, 2, 15075);
        tileBaseQualStats.addToCount(39, 0, 15063);
        tileBaseQualStats.addToCount(40, 0, 15073);
        tileBaseQualStats.addToCount(41, 2, 15063);
        tileBaseQualStats.addToCount(42, 3, 15043);
        tileBaseQualStats.addToCount(43, 0, 15089);
        tileBaseQualStats.addToCount(44, 2, 15088);
        tileBaseQualStats.addToCount(45, 1, 15058);
        tileBaseQualStats.addToCount(46, 2, 15078);
        tileBaseQualStats.addToCount(47, 0, 15034);
        tileBaseQualStats.addToCount(48, 1, 15038);
        tileBaseQualStats.addToCount(49, 2, 15011);
        tileBaseQualStats.addToCount(50, 2, 15034);
        tileBaseQualStats.addToCount(51, 1, 15017);
        tileBaseQualStats.addToCount(52, 2, 14995);
        tileBaseQualStats.addToCount(53, 4, 14979);
        tileBaseQualStats.addToCount(54, 2, 14964);
        tileBaseQualStats.addToCount(55, 3, 14937);
        tileBaseQualStats.addToCount(56, 1, 14971);
        tileBaseQualStats.addToCount(57, 1, 14968);
        tileBaseQualStats.addToCount(58, 0, 14942);
        tileBaseQualStats.addToCount(59, 4, 14962);
        tileBaseQualStats.addToCount(60, 3, 14934);
        tileBaseQualStats.addToCount(61, 566, 14303);
        tileBaseQualStats.addToCount(62, 1, 14914);
        tileBaseQualStats.addToCount(63, 2, 14954);
        tileBaseQualStats.addToCount(64, 1, 14922);
        tileBaseQualStats.addToCount(65, 2, 14875);
        tileBaseQualStats.addToCount(66, 1, 14896);
        tileBaseQualStats.addToCount(67, 1, 14847);
        tileBaseQualStats.addToCount(68, 0, 14811);
        tileBaseQualStats.addToCount(69, 0, 14810);
        tileBaseQualStats.addToCount(70, 2, 14892);
        tileBaseQualStats.addToCount(71, 0, 14799);
        tileBaseQualStats.addToCount(72, 2, 14831);
        tileBaseQualStats.addToCount(73, 1, 14772);
        tileBaseQualStats.addToCount(74, 1, 14829);
        tileBaseQualStats.addToCount(75, 0, 14775);
        tileBaseQualStats.addToCount(76, 2, 14882);
        tileBaseQualStats.addToCount(77, 1, 14825);
        tileBaseQualStats.addToCount(78, 2, 14813);
        tileBaseQualStats.addToCount(79, 508, 14323);
        tileBaseQualStats.addToCount(80, 1, 14687);
        tileBaseQualStats.addToCount(81, 0, 14729);
        tileBaseQualStats.addToCount(82, 0, 14664);
        tileBaseQualStats.addToCount(83, 2, 14785);
        tileBaseQualStats.addToCount(84, 0, 14729);
        tileBaseQualStats.addToCount(85, 4, 14659);
        tileBaseQualStats.addToCount(86, 2, 14691);
        tileBaseQualStats.addToCount(87, 1, 14677);
        tileBaseQualStats.addToCount(88, 2, 14641);
        tileBaseQualStats.addToCount(89, 2, 14520);
        tileBaseQualStats.addToCount(90, 0, 14621);
        tileBaseQualStats.addToCount(91, 1, 14629);
        tileBaseQualStats.addToCount(92, 2, 14457);
        tileBaseQualStats.addToCount(93, 1, 14588);
        tileBaseQualStats.addToCount(94, 1, 14596);
        tileBaseQualStats.addToCount(95, 0, 14523);
        tileBaseQualStats.addToCount(96, 1, 14579);
        tileBaseQualStats.addToCount(97, 1, 14576);
        tileBaseQualStats.addToCount(98, 1, 14565);
        tileBaseQualStats.addToCount(99, 2, 14572);
        tileBaseQualStats.addToCount(100, 0, 14549);
        tileBaseQualStats.addToCount(101, 2, 14484);
        tileBaseQualStats.addToCount(102, 1, 14554);
        tileBaseQualStats.addToCount(103, 3, 14456);
        tileBaseQualStats.addToCount(104, 1, 14426);
        tileBaseQualStats.addToCount(105, 1, 14400);
        tileBaseQualStats.addToCount(106, 3, 14477);
        tileBaseQualStats.addToCount(107, 2, 14397);
        tileBaseQualStats.addToCount(108, 2, 14332);
        tileBaseQualStats.addToCount(109, 0, 14347);
        tileBaseQualStats.addToCount(110, 0, 14322);
        tileBaseQualStats.addToCount(111, 0, 14266);
        tileBaseQualStats.addToCount(112, 1, 14255);
        tileBaseQualStats.addToCount(113, 2, 14331);
        tileBaseQualStats.addToCount(114, 3, 14105);
        tileBaseQualStats.addToCount(115, 3, 13996);
        tileBaseQualStats.addToCount(116, 4, 14198);
        tileBaseQualStats.addToCount(117, 0, 14220);
        tileBaseQualStats.addToCount(118, 0, 14052);
        tileBaseQualStats.addToCount(119, 1, 14113);
        tileBaseQualStats.addToCount(120, 3, 14144);
        tileBaseQualStats.addToCount(121, 2, 14159);
        tileBaseQualStats.addToCount(122, 1, 14087);
        tileBaseQualStats.addToCount(123, 3, 14077);
        tileBaseQualStats.addToCount(124, 4, 14099);
        tileBaseQualStats.addToCount(125, 1, 13993);
        tileBaseQualStats.addToCount(126, 1, 14038);
        tileBaseQualStats.addToCount(127, 0, 13994);
        tileBaseQualStats.addToCount(128, 3, 13999);
        tileBaseQualStats.addToCount(129, 2, 14024);
        tileBaseQualStats.addToCount(130, 1, 13975);
        tileBaseQualStats.addToCount(131, 3, 13946);
        tileBaseQualStats.addToCount(132, 3, 13876);
        tileBaseQualStats.addToCount(133, 5, 13785);
        tileBaseQualStats.addToCount(134, 1, 13715);
        tileBaseQualStats.addToCount(135, 0, 13757);
        tileBaseQualStats.addToCount(136, 0, 13757);
        tileBaseQualStats.addToCount(137, 2, 13716);
        tileBaseQualStats.addToCount(138, 0, 13660);
        tileBaseQualStats.addToCount(139, 0, 13591);
        tileBaseQualStats.addToCount(140, 1, 13439);
        tileBaseQualStats.addToCount(141, 2, 13492);
        tileBaseQualStats.addToCount(142, 2, 13548);
        tileBaseQualStats.addToCount(143, 1, 13411);
        tileBaseQualStats.addToCount(144, 2, 13309);
        tileBaseQualStats.addToCount(145, 1, 13439);
        tileBaseQualStats.addToCount(146, 1, 13505);
        tileBaseQualStats.addToCount(147, 0, 13505);
        tileBaseQualStats.addToCount(148, 2, 13491);
        tileBaseQualStats.addToCount(149, 3, 13405);
        tileBaseQualStats.addToCount(150, 10, 12654);



        /*
        TandemRepeat tandemRepeat = TandemRepeatFinder.findTandemRepeat("ATATATATCATAT");

        TestCase.assertEquals("AT", tandemRepeat.pattern);
        TestCase.assertEquals(8, tandemRepeat.length());

        tandemRepeat = TandemRepeatFinder.findTandemRepeat("ATGATGATGATCATATAT");

        TestCase.assertEquals("ATG", tandemRepeat.pattern);
        TestCase.assertEquals(11, tandemRepeat.length());

        tandemRepeat = TandemRepeatFinder.findTandemRepeat("TTTTTATGATGATGATCAGAGAGAGAGAGT");

        TestCase.assertEquals("AG", tandemRepeat.pattern);
        TestCase.assertEquals(12, tandemRepeat.length());

         */
    }
}
